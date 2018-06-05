//! Taxonomic tree index and "informativeness" lookup.

use cue::pipeline;
use daggy::Dag;
use daggy::petgraph::graph::NodeIndex;
use daggy::walker::Walker;
use error::*;

use index::TaxId;
use std::collections::{BTreeMap, BTreeSet};
use std::io::{BufRead, Write};
use std::process::exit;

/// Setting for determining how far "up" the tree to look for compatible siblings, based on the
/// number of parent edges to traverse.
#[derive(Clone, Copy, Debug)]
pub enum LcaSetting {
    /// Don't go up the tree at all -- this will only treat a result as informative if it only has
    /// one taxonomic ID matched.
    Zero,
    /// Go up one step in the tree (often to genus).
    One,
    /// Go up two steps in the tree (often to family).
    Two,
    /// Go up three steps in the tree.
    Three,
}

impl LcaSetting {
    /// Provide a list of numerical search heights to check in order.
    pub fn search_heights(&self) -> Vec<u32> {
        match *self {
            LcaSetting::Zero => vec![],
            LcaSetting::One => vec![1],
            LcaSetting::Two => vec![1, 2],
            LcaSetting::Three => vec![1, 2, 3],
        }
    }
}

/// Setting to determine how far "up" the tree to look for compatible siblings, based on the
/// logical rank of parent nodes.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum LogicalSetting {
    /// Collect all siblings that reside in the same genus.
    Genus,
    /// Collect all siblings that reside in the same genus or family, the latter occurring if the
    /// genus siblings are not informative.
    Family,
    /// Don't use logical ranks at all. This will force LCA querying for all reads.
    Na,
}

impl LogicalSetting {
    /// Provide a list of logical search heights to check in order.
    fn search_ranks(&self) -> Vec<Rank> {
        match *self {
            LogicalSetting::Genus => vec![Rank::Genus],
            LogicalSetting::Family => vec![Rank::Genus, Rank::Family],
            LogicalSetting::Na => vec![],
        }
    }
}

#[derive(RustcDecodable, RustcEncodable)]
struct Weight;

/// Index of the taxonomic tree.
///
/// Stores the tree as a DAG, and also contains a mapping from taxonomic IDs to DAG indices (the
/// DAG is stored in a flat array).
#[derive(RustcDecodable, RustcEncodable)]
pub struct TreeWithIndices {
    tree: Dag<TaxIdNode, Weight, u32>,
    indices: BTreeMap<TaxId, NodeIndex>,
}

impl TreeWithIndices {
    /// Taking an iterator of mtsv findings, determine which are "informative" according to the
    /// provided settings and write out a `read -> taxonomic ID` mapping. Executes in parallel.
    pub fn find_and_write_informatives<W: Write + Send + Sync>(&self,
                        all_hits: Box<Iterator<Item=MtsvResult<(String, BTreeSet<TaxId>)>>>,
                        num_threads: usize,
                        lca: LcaSetting,
                        logical: LogicalSetting,
                        output_wtr: &mut W) {

        pipeline("informative lookup pipeline",
                 num_threads,
                 all_hits,
                 |r| {

            // the iterator yields Result<...>
            if let Ok((read_id, hits)) = r {
                Some((read_id, self.informative_parent(&hits, lca, logical)))
            } else {
                error!("Problem parsing mtsv findings: {:?}", r);
                exit(1);
            }
        },
                 |result| {
            if let Some((read_id, Some(informative_tax_id))) = result {
                match write!(output_wtr, "{}:{}\n", read_id, informative_tax_id.0) {
                    Ok(_) => (),
                    Err(why) => {
                        error!("Problem writing informative results to file: {}", why);
                        exit(1);
                    },
                }
            }
        });
    }

    /// Determine an individual read's "informative parent": the lowest taxonomic node which has as
    /// children all taxonomic notes in the `hits` argument and which still conforms to the
    /// settings provided.
    fn informative_parent(&self,
                          hits: &BTreeSet<TaxId>,
                          lca: LcaSetting,
                          logical: LogicalSetting)
                          -> Option<TaxId> {
        // if there aren't any hits, nothing is informative
        if hits.len() == 0 {
            return None;
        }

        let mut ids = hits.iter();
        let mut first_id = *ids.next().unwrap();

        if hits.len() == 1 {
            // if there's only one hit, it's informative for that node by default,
            // no matter the settings
            return Some(first_id);
        }

        while let None = self.indices.get(&first_id) {
            first_id = match ids.next() {
                Some(id) => *id,
                None => {
                    error!("Unable to process a read -- all hits are missing from tax dump");
                    return None;
                },
            };
        }

        // search logical nodes first
        for rank in logical.search_ranks() {
            if let Some(sibs) = self.siblings_to_rank(first_id, rank) {
                let informed = self.get_informative(hits, &sibs);

                if informed.is_some() {
                    return informed;
                }
            }
        }

        // search LCA parents next
        for height in lca.search_heights() {
            let sibs = self.siblings_to_lca(first_id, height);
            let informed = self.get_informative(hits, &sibs);

            if informed.is_some() {
                return informed;
            }
        }

        None
    }

    fn get_informative(&self,
                       hits: &BTreeSet<TaxId>,
                       subindex: &(TaxId, BTreeSet<TaxId>))
                       -> Option<TaxId> {

        if hits.is_subset(&subindex.1) {
            Some(subindex.0)
        } else {
            None
        }
    }

    /// Find all "informative" siblings contained below the parent node at the specified logical
    /// rank.
    fn siblings_to_rank(&self, id: TaxId, rank: Rank) -> Option<(TaxId, BTreeSet<TaxId>)> {
        if let Some((genus, genus_idx)) = self.ancestor_at_rank(id, rank) {
            let mut rank_sibs = BTreeSet::new();
            rank_sibs.insert(genus);
            self.collect_children_to_depth(genus_idx, 100, &mut rank_sibs);

            Some((genus, rank_sibs))
        } else {
            None
        }
    }

    /// Find all "informative" siblings contained below the parent node at the specified "edge
    /// height" (number of steps up).
    fn siblings_to_lca(&self, id: TaxId, lca: u32) -> (TaxId, BTreeSet<TaxId>) {
        let (parent, parent_idx) = self.ancestor_at_height(id, lca);
        let mut lca_sibs = BTreeSet::new();
        lca_sibs.insert(parent);
        self.collect_children_to_depth(parent_idx, 100, &mut lca_sibs);

        (parent, lca_sibs)
    }

    /// Find the parent of a node at a given logical rank, limited to 25 jumps in cases where a
    /// node has no parent at that rank.
    fn ancestor_at_rank(&self, id: TaxId, rank: Rank) -> Option<(TaxId, NodeIndex)> {
        let idx = match self.indices.get(&id) {
            Some(idx) => *idx,
            None => {
                error!("Asking for tax ID {:?} which is not present in tax dump.",
                       id);
                exit(1);
            },
        };

        self.rank_ancestor_helper(idx, rank, 30)
    }

    /// Recursively move up the tree, looking for a parent at a given Rank.
    fn rank_ancestor_helper(&self,
                            idx: NodeIndex,
                            rank: Rank,
                            remain: u32)
                            -> Option<(TaxId, NodeIndex)> {
        if remain == 0 {
            return None;
        }

        self.tree.node_weight(idx).and_then(|w| {
            if w.rank == rank {
                Some((w.id, idx))
            } else {
                self.tree
                    .parents(idx)
                    .next_node(&self.tree)
                    .and_then(|p| self.rank_ancestor_helper(p, rank, remain - 1))
            }
        })
    }

    /// Find the parent node at a given "height" (number of edge jumps).
    fn ancestor_at_height(&self, id: TaxId, height: u32) -> (TaxId, NodeIndex) {

        let child_index = match self.indices.get(&id) {
            Some(idx) => *idx,
            None => {
                error!("Asking for a tax ID ({:?}) which is not present in the tax dump.",
                       id);
                exit(1);
            },
        };

        // we'll naively assume that all interesting leaf nodes have parents up to the
        // specified height
        let mut parent = child_index;
        let mut parent_id = id;
        for _ in 0..height {

            match self.tree.parents(parent).next_node(&self.tree) {
                Some(p) => parent = p,
                None => break,
            }

            parent_id = self.tree.node_weight(parent).unwrap().id;
        }

        (parent_id, parent)
    }

    /// Collect all children of a given tree node, up to a certain depth.
    fn collect_children_to_depth(&self,
                                 idx: NodeIndex,
                                 height: u32,
                                 to_add_to: &mut BTreeSet<TaxId>) {
        // if there are in fact children but we've reached our maximum depth
        if height == 0 {
            return;
        }

        // if there are no children, this will gracefully noop
        for (_, child_index) in self.tree.children(idx).iter(&self.tree) {

            let child_taxid = self.tree.node_weight(child_index).unwrap().id;
            to_add_to.insert(child_taxid);

            self.collect_children_to_depth(child_index, height - 1, to_add_to);
        }
    }

    /// Construct a taxonomic tree index from a filename.
    pub fn from_node_dump<R: BufRead>(node_dump: R) -> MtsvResult<Self> {
        info!("Found nodes.dmp, tokenizing and lexing taxon nodes...");
        let (children, ranks) = try!(Self::tokenize_nodes(node_dump));

        info!("Finished lexing nodes, parsing taxonomic tree...");
        Ok(TreeWithIndices::from_loose_nodes(children, ranks))
    }

    /// Construct a taxonomic tree index from a series of nodes with specified children but no
    /// structure.
    fn from_loose_nodes(children: BTreeMap<TaxId, Vec<TaxId>>,
                        ranks: BTreeMap<TaxId, Rank>)
                        -> Self {

        let mut tree = Dag::<TaxIdNode, Weight, u32>::new();
        let mut indices = BTreeMap::new();

        let default_rank = Rank::Other;

        for (parent, children) in children.into_iter() {

            // if the parent's already insterted from another parent's child list
            let parent_idx = indices.get(&parent).map(|m| *m).unwrap_or_else(|| {
                // otherwise this node hasn't been inserted yet
                let rank = ranks.get(&parent).unwrap_or(&default_rank);
                let new_node = TaxIdNode {
                    id: parent,
                    rank: *rank,
                };
                let idx = tree.add_node(new_node);
                indices.insert(parent, idx);
                idx
            });

            // now the current parent is in the graph, and we have its index
            for child_taxid in children {
                if !indices.contains_key(&child_taxid) {
                    let rank = ranks.get(&child_taxid).unwrap_or(&default_rank);

                    let new_node = TaxIdNode {
                        id: child_taxid,
                        rank: *rank,
                    };

                    indices.insert(child_taxid, tree.add_child(parent_idx, Weight, new_node).1);
                } else {
                    // unwrap is OK, b/c the map already has the child here
                    let _ = tree.add_edge(parent_idx, *indices.get(&child_taxid).unwrap(), Weight);
                }
            }
        }

        TreeWithIndices {
            tree: tree,
            indices: indices,
        }
    }

    /// Parse the node dump file from NCBI into structures of loose nodes.
    fn tokenize_nodes<R: BufRead>
        (contents: R)
         -> MtsvResult<(BTreeMap<TaxId, Vec<TaxId>>, BTreeMap<TaxId, Rank>)> {

        let mut ranks = BTreeMap::new();
        let mut children = BTreeMap::new();

        for line in contents.lines() {
            let line = try!(line);
            let trimmed = line.trim();

            if trimmed.len() == 0 {
                continue; // we shouldn't care about pure whitespace lines
            }

            let mut tokens = trimmed.split("|").map(|t| t.trim());

            let taxid = match tokens.next().map(|t| t.parse::<TaxId>()) {
                Some(Ok(t)) => t,
                _ => return Err(MtsvError::InvalidHeader(trimmed.to_string())),
            };

            let parentid = match tokens.next().map(|t| t.parse::<TaxId>()) {
                Some(Ok(t)) => t,
                _ => return Err(MtsvError::InvalidHeader(trimmed.to_string())),
            };

            let rank = match tokens.next() {
                Some(t) => Rank::from(t),
                None => return Err(MtsvError::InvalidHeader(trimmed.to_string())),
            };

            children.entry(parentid).or_insert(Vec::new()).push(taxid);
            ranks.insert(taxid, rank);
        }

        Ok((children, ranks))
    }
}

#[allow(missing_docs)]
#[derive(Copy, Clone, RustcDecodable, RustcEncodable)]
pub struct TaxIdNode {
    id: TaxId,
    rank: Rank,
}

#[allow(missing_docs)]
#[derive(Copy, Clone, Debug, Eq, PartialEq, RustcDecodable, RustcEncodable)]
pub enum Rank {
    SpeciesSubgroup,
    SubFamily,
    Phylum,
    SubClass,
    SubPhylum,
    Family,
    SubOrder,
    Kingdom,
    InfraOrder,
    Order,
    SuperPhylum,
    Forma,
    SuperKingdom,
    SpeciesGroup,
    Varietas,
    ParvOrder,
    InfraClass,
    SubSpecies,
    SuperFamily,
    SuperClass,
    Tribe,
    SubTribe,
    SubGenus,
    Class,
    Species,
    Genus,
    SubKingdom,
    SuperOrder,
    NoRank,
    Other,
}

impl<'a> From<&'a str> for Rank {
    fn from(s: &'a str) -> Self {
        match s {
            "species subgroup" => Rank::SpeciesSubgroup,
            "subfamily" => Rank::SubFamily,
            "phylum" => Rank::Phylum,
            "subclass" => Rank::SubClass,
            "subphylum" => Rank::SubPhylum,
            "family" => Rank::Family,
            "suborder" => Rank::SubOrder,
            "kingdom" => Rank::Kingdom,
            "infraorder" => Rank::InfraOrder,
            "order" => Rank::Order,
            "superphylum" => Rank::SuperPhylum,
            "forma" => Rank::Forma,
            "superkingdom" => Rank::SuperKingdom,
            "species group" => Rank::SpeciesGroup,
            "varietas" => Rank::Varietas,
            "parvorder" => Rank::ParvOrder,
            "infraclass" => Rank::InfraClass,
            "subspecies" => Rank::SubSpecies,
            "superfamily" => Rank::SuperFamily,
            "superclass" => Rank::SuperClass,
            "tribe" => Rank::Tribe,
            "subtribe" => Rank::SubTribe,
            "subgenus" => Rank::SubGenus,
            "class" => Rank::Class,
            "species" => Rank::Species,
            "genus" => Rank::Genus,
            "subkingdom" => Rank::SubKingdom,
            "superorder" => Rank::SuperOrder,
            "no rank" => Rank::NoRank,
            _ => Rank::Other,
        }
    }
}

#[cfg(test)]
mod test {
    use index::TaxId;
    use std::collections::{BTreeMap, BTreeSet};
    use std::io::Cursor;
    use super::*;

    #[test]
    fn informative_parent() {
        let tree = construct();

        let mut findings = BTreeSet::new();

        assert_eq!(None,
                   tree.informative_parent(&findings, LcaSetting::Three, LogicalSetting::Family));

        findings.insert(TaxId(9606));

        assert_eq!(Some(TaxId(9606)),
                   tree.informative_parent(&findings, LcaSetting::Zero, LogicalSetting::Na));

        findings.insert(TaxId(1425170));

        assert_eq!(None,
                   tree.informative_parent(&findings, LcaSetting::Zero, LogicalSetting::Na));

        assert_eq!(Some(TaxId(9605)),
                   tree.informative_parent(&findings, LcaSetting::One, LogicalSetting::Na));

        assert_eq!(Some(TaxId(9605)),
                   tree.informative_parent(&findings, LcaSetting::Zero, LogicalSetting::Genus));

        assert_eq!(Some(TaxId(9605)),
                   tree.informative_parent(&findings, LcaSetting::One, LogicalSetting::Genus));

        findings.insert(TaxId(9592));

        assert_eq!(None,
                   tree.informative_parent(&findings, LcaSetting::Zero, LogicalSetting::Na));

        assert_eq!(None,
                   tree.informative_parent(&findings, LcaSetting::One, LogicalSetting::Na));

        assert_eq!(None,
                   tree.informative_parent(&findings, LcaSetting::Zero, LogicalSetting::Genus));

        assert_eq!(None,
                   tree.informative_parent(&findings, LcaSetting::One, LogicalSetting::Genus));

        assert_eq!(Some(TaxId(9604)),
                   tree.informative_parent(&findings, LcaSetting::Two, LogicalSetting::Na));

        assert_eq!(Some(TaxId(9604)),
                   tree.informative_parent(&findings, LcaSetting::Zero, LogicalSetting::Family));

        assert_eq!(Some(TaxId(9604)),
                   tree.informative_parent(&findings, LcaSetting::Two, LogicalSetting::Family));

        findings.insert(TaxId(29089));

        assert_eq!(None,
                   tree.informative_parent(&findings, LcaSetting::Zero, LogicalSetting::Na));

        assert_eq!(None,
                   tree.informative_parent(&findings, LcaSetting::One, LogicalSetting::Na));

        assert_eq!(None,
                   tree.informative_parent(&findings, LcaSetting::Zero, LogicalSetting::Genus));

        assert_eq!(None,
                   tree.informative_parent(&findings, LcaSetting::One, LogicalSetting::Genus));

        assert_eq!(None,
                   tree.informative_parent(&findings, LcaSetting::Two, LogicalSetting::Na));

        assert_eq!(None,
                   tree.informative_parent(&findings, LcaSetting::Zero, LogicalSetting::Family));

        assert_eq!(None,
                   tree.informative_parent(&findings, LcaSetting::Two, LogicalSetting::Family));

        assert_eq!(Some(TaxId(314295)),
                   tree.informative_parent(&findings, LcaSetting::Three, LogicalSetting::Na));
    }

    #[test]
    fn human_siblings_rank() {
        let tree = construct();

        let mut siblings = BTreeSet::new();
        let leaf = TaxId(716692);

        let (found_inform, found_sibs) = tree.siblings_to_rank(leaf, Rank::SubSpecies).unwrap();
        siblings.insert(leaf);
        assert_eq!(found_inform, leaf);
        assert_eq!(&found_sibs, &siblings);

        let (found_inform, found_sibs) = tree.siblings_to_rank(leaf, Rank::Species).unwrap();
        siblings.extend(vec![TaxId(9588), TaxId(716694), TaxId(716693)]);
        assert_eq!(found_inform, TaxId(9588));
        assert_eq!(&found_sibs, &siblings);

        let (found_inform, found_sibs) = tree.siblings_to_rank(leaf, Rank::Genus).unwrap();
        siblings.extend(vec![TaxId(9578),
                             TaxId(288892),
                             TaxId(9581),
                             TaxId(9589),
                             TaxId(81572),
                             TaxId(9587),
                             TaxId(716691),
                             TaxId(9579),
                             TaxId(716690),
                             TaxId(9580),
                             TaxId(716695),
                             TaxId(716696),
                             TaxId(716697),
                             TaxId(716698)]);
        assert_eq!(found_inform, TaxId(9578));
        assert_eq!(&found_sibs, &siblings);

        let (found_inform, found_sibs) = tree.siblings_to_rank(leaf, Rank::Family).unwrap();
        siblings.extend(vec![TaxId(9577),
                             TaxId(325165),
                             TaxId(1616038),
                             TaxId(61852),
                             TaxId(693984),
                             TaxId(61853),
                             TaxId(327374),
                             TaxId(9586),
                             TaxId(29089),
                             TaxId(716684),
                             TaxId(423450),
                             TaxId(423449),
                             TaxId(716685),
                             TaxId(325166),
                             TaxId(9590),
                             TaxId(405039),
                             TaxId(325167),
                             TaxId(61851),
                             TaxId(593543)]);
        assert_eq!(found_inform, TaxId(9577));
        assert_eq!(&found_sibs, &siblings);
    }

    #[test]
    fn human_siblings_lca() {
        let tree = construct();

        let mut siblings = BTreeSet::new();
        let leaf = TaxId(46359);

        let (found_inform, found_sibs) = tree.siblings_to_lca(leaf, 0);
        siblings.insert(leaf);
        assert_eq!(found_inform, leaf);
        assert_eq!(&found_sibs, &siblings);

        let (found_inform, found_sibs) = tree.siblings_to_lca(leaf, 1);
        siblings.insert(TaxId(1159185));
        assert_eq!(found_inform, TaxId(1159185));
        assert_eq!(&found_sibs, &siblings);

        let (found_inform, found_sibs) = tree.siblings_to_lca(leaf, 2);
        siblings.insert(TaxId(499232));
        assert_eq!(found_inform, TaxId(499232));
        assert_eq!(&found_sibs, &siblings);

        let new_lca3_sibs =
            vec![TaxId(9592), TaxId(9593), TaxId(182511), TaxId(9595), TaxId(406788)];
        siblings.extend(new_lca3_sibs);
        let (found_inform, found_sibs) = tree.siblings_to_lca(leaf, 3);
        assert_eq!(found_inform, TaxId(9592));
        assert_eq!(&found_sibs, &siblings);

        siblings.insert(TaxId(207598));
        let (found_inform, found_sibs) = tree.siblings_to_lca(leaf, 4);
        assert_eq!(found_inform, TaxId(207598));
        assert_eq!(&found_sibs, &siblings);

        let new_lca5_sibs = vec![TaxId(9604),
                                 TaxId(9596),
                                 TaxId(9597),
                                 TaxId(9598),
                                 TaxId(756884),
                                 TaxId(37010),
                                 TaxId(37011),
                                 TaxId(91950),
                                 TaxId(37012),
                                 TaxId(1294088),
                                 TaxId(9605),
                                 TaxId(1425170),
                                 TaxId(9606),
                                 TaxId(63221),
                                 TaxId(741158)];
        siblings.extend(new_lca5_sibs);
        let (found_inform, found_sibs) = tree.siblings_to_lca(leaf, 5);
        assert_eq!(found_inform, TaxId(9604));
        assert_eq!(&found_sibs, &siblings);
    }

    #[test]
    fn human_ancestors() {
        // see tests/nodes.dmp for the structure here
        let tree = construct();

        assert_eq!(TaxId(182511),
                   tree.ancestor_at_rank(TaxId(182511), Rank::SubSpecies).unwrap().0);
        assert_eq!(TaxId(9593),
                   tree.ancestor_at_rank(TaxId(182511), Rank::Species).unwrap().0);
        assert_eq!(TaxId(9592),
                   tree.ancestor_at_rank(TaxId(182511), Rank::Genus).unwrap().0);
        assert_eq!(TaxId(207598),
                   tree.ancestor_at_rank(TaxId(182511), Rank::SubFamily).unwrap().0);
        assert_eq!(TaxId(9604),
                   tree.ancestor_at_rank(TaxId(182511), Rank::Family).unwrap().0);
        assert_eq!(TaxId(314295),
                   tree.ancestor_at_rank(TaxId(182511), Rank::SuperFamily).unwrap().0);

        // the root node will always be other, but it should be the only one
        assert_eq!(TaxId(1),
                   tree.ancestor_at_rank(TaxId(182511), Rank::Other).unwrap().0);

        // make sure we get none out
        assert_eq!(None,
                   tree.ancestor_at_rank(TaxId(182511), Rank::SpeciesSubgroup));
        assert_eq!(None, tree.ancestor_at_rank(TaxId(182511), Rank::Varietas));
        assert_eq!(None, tree.ancestor_at_rank(TaxId(182511), Rank::Tribe));
        assert_eq!(None, tree.ancestor_at_rank(TaxId(182511), Rank::SubTribe));

        assert_eq!(TaxId(9604), tree.ancestor_at_height(TaxId(9604), 0).0);
        assert_eq!(TaxId(314295), tree.ancestor_at_height(TaxId(9604), 1).0);
        assert_eq!(TaxId(9526), tree.ancestor_at_height(TaxId(9604), 2).0);
        assert_eq!(TaxId(314293), tree.ancestor_at_height(TaxId(9604), 3).0);
        assert_eq!(TaxId(376913), tree.ancestor_at_height(TaxId(9604), 4).0);
        assert_eq!(TaxId(9443), tree.ancestor_at_height(TaxId(9604), 5).0);
        assert_eq!(TaxId(314146), tree.ancestor_at_height(TaxId(9604), 6).0);

        assert_eq!(TaxId(406788), tree.ancestor_at_height(TaxId(406788), 0).0);
        assert_eq!(TaxId(9593), tree.ancestor_at_height(TaxId(406788), 1).0);
        assert_eq!(TaxId(9592), tree.ancestor_at_height(TaxId(406788), 2).0);
        assert_eq!(TaxId(207598), tree.ancestor_at_height(TaxId(406788), 3).0);
        assert_eq!(TaxId(9604), tree.ancestor_at_height(TaxId(406788), 4).0);
        assert_eq!(TaxId(314295), tree.ancestor_at_height(TaxId(406788), 5).0);

        assert_eq!(TaxId(1), tree.ancestor_at_height(TaxId(406788), 1_000).0);
    }

    fn construct() -> TreeWithIndices {
        let to_parse = include_str!("../tests/nodes.dmp");
        TreeWithIndices::from_node_dump(Cursor::new(to_parse)).unwrap()
    }

    #[test]
    fn tokenize_empty() {
        let found = TreeWithIndices::tokenize_nodes(Cursor::new(String::from(""))).unwrap();

        assert_eq!(found, (BTreeMap::new(), BTreeMap::new()));
    }

    #[test]
    #[should_panic]
    fn tokenize_fail_taxid_bad() {
        let line = "abc | 122345 | species".to_string();
        TreeWithIndices::tokenize_nodes(Cursor::new(line)).unwrap();
    }

    #[test]
    #[should_panic]
    fn tokenize_fail_taxid_missing() {
        let line = " | 122345 | species".to_string();
        TreeWithIndices::tokenize_nodes(Cursor::new(line)).unwrap();
    }

    #[test]
    #[should_panic]
    fn tokenize_fail_parent_bad() {
        let line = "12345 | abc | species".to_string();
        TreeWithIndices::tokenize_nodes(Cursor::new(line)).unwrap();
    }

    #[test]
    #[should_panic]
    fn tokenize_fail_parent_missing() {
        let line = "12345 |  | species".to_string();
        TreeWithIndices::tokenize_nodes(Cursor::new(line)).unwrap();
    }

    #[test]
    #[should_panic]
    fn tokenize_fail_rank_missing() {
        let line = "1234567 | 122345".to_string();
        TreeWithIndices::tokenize_nodes(Cursor::new(line)).unwrap();
    }

    #[test]
    fn tokenize() {
        let lines = "9606    |       9605    |       species |       HS      |       5      \
         |       1       |       1       |       1       |       2  |1       |       1      \
          |       0       |               |
96060   |       141711  |       genus   |               |       1       |       1       |      \
 1       |       1       |       5  |1       |       0       |       0       |               |
 1425170 |       9605    |       species |       HH      |       5       |       1       |      \
  1       |       1       |       2  |1       |       1       |       0       |               |\
";

        let mut ranks = BTreeMap::new();
        ranks.insert(TaxId(9606), Rank::Species);
        ranks.insert(TaxId(1425170), Rank::Species);
        ranks.insert(TaxId(96060), Rank::Genus);

        let mut parents = BTreeMap::new();
        parents.insert(TaxId(9605), vec![TaxId(9606), TaxId(1425170)]);
        parents.insert(TaxId(141711), vec![TaxId(96060)]);

        let found = TreeWithIndices::tokenize_nodes(Cursor::new(lines)).unwrap();

        assert_eq!(found, (parents, ranks));
    }
}
