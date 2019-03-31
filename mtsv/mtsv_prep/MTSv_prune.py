import subprocess, shutil
from io import BytesIO
import argparse
import fnmatch
import os, datetime
from ftplib import FTP, all_errors
from time import sleep
import gzip
import tarfile
from multiprocessing import Pool, Queue, Process, Manager, RLock
import pickle, json

lock = RLock()

def serialization(gi2tx,fasta_path, txdump_path):

    tx2gi = {}
    gi2ind = {}
    print("Parsing taxid to Unique ID")

    for i in gi2tx:
        try:
            with gzip.open(i) as file:
                for line in file:
                    line = line.strip().split(b'\t')
                    try:
                        tx2gi[line[1].strip()].append(line[0].strip())
                    except:
                        tx2gi[line[1].strip()] = [line[0].strip()]
        except:
            with open(i, "rb") as file:
                for line in file:
                    line = line.strip().split(b'\t')
                    try:
                        tx2gi[line[1].strip()].append(line[0].strip())
                    except:
                        tx2gi[line[1].strip()] = [line[0].strip()]
        print("Indexing Fasta")
        with open(fasta_path, "rb") as file:
            for line in file:
                if chr(line[0]) == ">":
                    gi = line.split(b' ',2)[1].split(b':')[1].strip()
                    gi2ind[gi] = file.tell()-len(line)
    print("Mapping taxid to index")
    for key in tx2gi.keys():
        temp = []
        for gi in tx2gi[key]:
            try:
                 temp.append(gi2ind[gi])
            except KeyError:
                continue
        tx2gi[key] = temp
    tx_ids, child2parent = taxids2name(txdump_path)

    print("Serializing")
    out = os.path.abspath(fasta_path).rsplit(".",1)[0]+".p"
    with open(out, "wb") as file:
        pickle.dump([tx_ids, child2parent,tx2gi], file)

def acc_serialization(acc2tx,fasta_path, txdump_path, overwrite=True):

    tx2gi = {}
    acc2ind = {}
    # print("Indexing Fasta")
    with open(fasta_path, "rb") as file:
        for line in file:
            if chr(line[0]) == ">":
                acc = line.split(b' ', 1)[0][1:].split(b'.')[0].strip()
                acc2ind[acc] = file.tell() - len(line)


    # print("Parsing taxid to Unique ID")

    for i in acc2tx:
        try:
            with gzip.open(i) as file:
                for line in file:
                    line = line.strip().split(b'\t')
                    try:
                        acc2ind[line[0].strip()]
                    except KeyError:
                        continue
                    try:
                        tx2gi[line[2].strip()].append(line[0].strip())
                    except:
                        tx2gi[line[2].strip()] = [line[0].strip()]
        except:
            with open(i, "rb") as file:
                for line in file:
                    line = line.strip().split(b'\t')
                    try:
                        acc2ind[line[0].strip()]
                    except KeyError:
                        continue
                    try:
                        tx2gi[line[2].strip()].append(line[0].strip())
                    except:
                        tx2gi[line[2].strip()] = [line[0].strip()]
    # print("Mapping taxid to index")
    for key in tx2gi.keys():
        temp = []
        for gi in tx2gi[key]:
            try:
                 temp.append(acc2ind[gi])
            except KeyError:
                continue
        tx2gi[key] = temp
    tx_ids, child2parent = taxids2name(txdump_path)

    # print("Serializing")
    out = os.path.abspath(fasta_path).rsplit(".",1)[0]+".p"
    with open(out, "wb") as file:
        pickle.dump([tx_ids, child2parent,tx2gi], file)
    return out

def deserialization(pickle_path):
    with open(pickle_path, "rb") as file:
        return pickle.load(file)

# Perform a binary search of a list, returning the found index or -1 if not
def binary_search(a, x, lo=0, hi=None):   # can't use a to specify default for hi
    hi = len(a)-1
    lo = 0
    mid = (hi-lo)//2
    while True:
        if lo >= hi:
            if a[mid][0] == x: return mid
            else: return -1
        if a[mid][0] == x:
            return mid
        elif a[mid][0] < x:
            lo = mid +1
        else:
            hi = mid -1
        mid = ((hi - lo) // 2)+lo

def taxids2name(dump_path):
    dump = tarfile.open(dump_path, "r:gz")
    tax_ids = {}
    names = dump.extractfile(dump.getmember('nodes.dmp')).read().split(b'\t|\n')
    roller = {}
    for l in names:
        if len(l) == 0:
            continue
        tokens = l.split(b'\t|\t')

        c_taxid = tokens[0].decode().strip()
        p_taxid = tokens[1].decode().strip()
        roller[tokens[0]] = (tokens[1],tokens[2])
        try:
            tax_ids[c_taxid]
        except KeyError:
            tax_ids[c_taxid] = []
        try:
            if p_taxid != c_taxid:
                tax_ids[p_taxid].append(c_taxid)
        except:
            tax_ids[p_taxid] = [c_taxid]

    return tax_ids, roller

# Performs breadth first search of NCBI tree returning a set of leaf tax nodes with offsets in the fasta DB
def get_tree(tx_ids, tx, terminals):
    # try:
    #     terminals[tx.encode()]
    nodes = {tx}
    # except KeyError:
    # nodes = set()
    try:
        cur_level = tx_ids[tx]
    except KeyError:
        return nodes
    temp = []

    while True:
        for tx in cur_level:
            try:
                terminals[tx.encode()]
                nodes.add(tx)
            except KeyError:
                pass
            finally:
                temp += tx_ids[tx]
        if not temp:
            break
        cur_level = temp
        temp = []

    return nodes

# This is used to roll a NCBI taxonomic ID to a desired rank
def roll_up(tx_id, rank, c2p, prev_roll=None):
    if tx_id:
        try:
            if prev_roll:
                return prev_roll[tx_id]
            else:
                raise KeyError
        except KeyError:
            try:
                cur = tx_id
                # nxt = c2p[cur][0]
                while cur != b'1':
                    # print(cur)
                    if c2p[cur][1] == rank:
                        if prev_roll:
                            prev_roll[tx_id] = cur
                        # print(tx_id)
                        return cur
                    else:
                        cur = c2p[cur][0]
                if tx_id and cur == b'1':
                    prev_roll[tx_id] = tx_id
            except TypeError:
                pass
        return tx_id

        # else:
        #     return
# This is the function that will extract the taxonomic sequences of interest from the parsed flat file fasta
# process is:
# Parse inclusive and exclusive taxids
# deserialize previously built NCBI tree structure and byte offset information
# Call function (uses depth first search) to get a set of all child of the taxid repeat for tax ids to exclude
# use set difference to get desired leaf taxids
# Opens fasta DB and out file reading sequence in from start of sequence header roll up occuring at runtime
def clip(in_tx,ru_rank, ex_tx, name, min,maximum,fasta_path, pickle_path, debug=False):
    if debug:
        return os.path.abspath(name)
    # try:
    if len(in_tx) ==1:
        temp = []
        try:
            with open(in_tx[0], "r") as file:
                for line in file:
                    temp.append(line.strip())
            in_tx = temp
        except FileNotFoundError:
            pass
    if ex_tx and len(ex_tx) ==1:
        temp = []
        try:
            with open(ex_tx[0],"r") as file:
                for line in file:
                    temp.append(line.strip())
            ex_tx = temp
        except FileNotFoundError:
            pass
    # print(pickle_path)
    # print("Getting Offsets and Tree")
    tx_ids, child2parent, positions = deserialization(pickle_path)

    taxons = set()
    # print("Getting TaxIds")
    for i in in_tx:
        taxons = taxons.union(get_tree(tx_ids, i, positions))

    if ex_tx:
        # print("\tPruning")
        for i in ex_tx:
           taxons = taxons.difference(get_tree(tx_ids, i, positions))

    if not name:
        name = "_".join(in_tx)
        if ex_tx:
            name += "_not_{0}".format("_".join(ex_tx))
        name += "_seqs.fasta"
    if ru_rank:
        ru_rank = ru_rank.encode()

    seq = bytearray()
    line_count = 0
    # print("Writing")
    with open(fasta_path, "rb") as fasta:
        with open(name, "wb") as out:
            for tx in sorted(taxons):
                if tx:
                    tx = tx.encode().strip()
                    try:
                        positions[tx].sort()
                    except KeyError:
                        continue
                    if ru_rank:
                        rr_tx = roll_up(tx, ru_rank, child2parent)
                    else:
                        rr_tx = tx
                    if not rr_tx:
                        continue
                    # print()
                    positions[tx].sort()
                    for off in positions[tx]:
                        fasta.seek(off)
                        header = fasta.readline()
                        line = fasta.readline()
                        while line and chr(line[0]) != ">":
                            seq += line
                            line_count += 1
                            line = fasta.readline()
                        if len(seq)-line_count >= min and len(seq)-float(line_count)<= maximum:
                            gi = header.split(b' ',2)[1].split(b':')[1]
                            # gi = header.split(b' ',1)[0].strip(b'>')
                            out.write(b'>'+gi+b'-'+rr_tx+b'\n')
                            out.write(seq)
                        line_count = 0
                        seq = bytearray()
    return os.path.abspath(name)
# except:
    #     return 0
# writes a new or updates json config file
def gen_json(configuration, args):
    if args.update and args.configuration_path:
        with open(args.configuration_path, "w") as file:
            json.dump(configuration, file, sort_keys=True)

    else:
        with open(args.output + ".json", "w") as file:
            json.dump(configuration, file, sort_keys=True, indent=4)

# Returns a dictionary of the config file
def parse_json(args):
    try:
        with open(args.configuration_path, "r") as file:
            return json.load(file)
    except:
        with open(args, "r") as file:
            return json.load(file)

# Parses command line arguments into dictionary of arguments or adds to one.  Can be serialized
def arg_unwrappers(args, arguments=None):
    if not arguments:
        arguments = {}

    if args.serialization_path:
        arguments['serialization-path'] = os.path.abspath( args.serialization_path)
    try:
        arguments['serialization-path']
    except KeyError:
        arguments['serialization-path'] = ""

    if args.fasta_path:
        arguments['fasta-path'] = os.path.abspath(args.fasta_path)
    try:
        arguments['fasta-path']
    except KeyError:
        arguments['fasta-path'] = ""

    if args.minimum_length:
        arguments['minimum-length'] = args.minimum_length
    try:
        arguments['minimum-length']
    except KeyError:
        arguments['minimum-length'] = 0

    if args.maximum_length:
        arguments['maximum-length'] = args.maximum_length
    try:
        arguments['maximum-length']
    except KeyError:
        arguments['maximum-length'] = float('inf')
    if args.taxdump_path:
        arguments['taxdump-path'] = os.path.abspath(args.taxdump_path)
    try:
        arguments['taxdump-path']
    except:
        arguments['taxdump-path'] = ""

    # if args.gi_to_taxid_paths:
    #     temp = []
    #     for i in args.gi_to_taxid_path:
    #         temp.append(os.path.abspath(i))
    #
    #     arguments['gi-to-taxid-paths'] = temp
    # try:
    #     arguments['gi-to-taxid-paths']
    # except:
    #     arguments['gi-to-taxid-paths'] = []

    if args.acc_to_taxid_paths:
        temp = []
        for i in args.acc_to_taxid_paths:
            temp.append(os.path.abspath(i))

        arguments['acc-to-taxid-paths'] = temp
    try:
        arguments['acc-to-taxid-paths']
    except:
        arguments['acc-to-taxid-paths'] = []

    if args.rollup_rank:
        arguments['rollup-rank'] = args.rollup_rank.lower()
    try:
        arguments['rollup-rank']
    except:
        arguments['rollup-rank'] = ""


    return arguments

def oneclickjson(path):
    arguments = []
    for fh in os.listdir(os.path.join(path, "artifacts")):
        if fnmatch.fnmatch(fh, "*.fas"):
            arguments.append({})
            fh = fh.split(".")[0]
            arguments[-1]['serialization-path'] = os.path.abspath( os.path.join(path, "artifacts","{0}.p".format(fh)))

            arguments[-1]['fasta-path'] = os.path.abspath(os.path.join(path, "artifacts", "{0}.fas".format(fh)))

            arguments[-1]['minimum-length'] = 0

            arguments[-1]['maximum-length'] = float('inf')
            arguments[-1]['taxdump-path'] = os.path.abspath(os.path.join(path, "artifacts","taxdump.tar.gz"))

            arguments[-1]['acc-to-taxid-paths'] = []
            for fp in os.listdir(os.path.join(path,"artifacts")):
                    if fnmatch.fnmatch(fp, "*accession2taxid*"):
                        arguments[-1]['acc-to-taxid-paths'].append(os.path.abspath(os.path.join(path,"artifacts",fp)))

            arguments[-1]['rollup-rank'] = "species"

            with open(os.path.abspath(os.path.join(path,"artifacts", "{0}.json".format(fh))), "w") as file:
                json.dump(arguments[-1], file, sort_keys=True, indent=4)

    return arguments

def build_db( flat_list_in_fp, fasta_out_fp, keyword_out_fp, source_out_fp, thread_count, gi_to_word):
    start_dir = os.getcwd()
    script_dir = os.path.dirname(os.path.realpath(__file__))
    flat_list_in_fp = os.path.abspath(flat_list_in_fp)

    os.chdir(script_dir)
    command_one = "g++ -std=c++11 -pthread -static-libstdc++ taxidtool.cpp -o db_builder"
    command_two = "./db_builder {0} {1}.tmp {2}.tsv {3}".format(flat_list_in_fp, fasta_out_fp, fasta_out_fp.rsplit(".",1)[0],
                                                                        thread_count)

    command_three = "mtsv-db-build {0} {1}.tmp {2}.tsv {3}".format(flat_list_in_fp, fasta_out_fp, fasta_out_fp.rsplit(".",1)[0],
                                                                     thread_count)
    if not os.path.isfile(fasta_out_fp+".tmp") and not os.path.isfile(fasta_out_fp):
        try:
            subprocess.run(command_three.split())
        except:
            try:
                subprocess.run(command_two.split())
            except:
                subprocess.run(command_one.split())
                subprocess.run(command_two.split())

    count = 0
    if os.path.isfile(fasta_out_fp+".tmp") and not os.path.isfile(fasta_out_fp):
        with open("{0}.tmp".format(fasta_out_fp), "rb") as start_file:
            with open(fasta_out_fp+".temp", "wb") as end_file:
                for line in start_file:
                    if chr(line[0]) == ">":
                        line = line.split(b'GI:')
                        if len(line) == 2:
                            line[1] = line[1].split(b' ',1)[1].strip()+b'\n'
                            line[0] = line[0].strip()
                            header = line
                        else:
                            header = line[0].split(b' ', 1)
                        end_file.write(" GI:{0} ".format(count).encode().join(header))
                        count += 1
                    else:
                        end_file.write(line)
        os.rename(fasta_out_fp+".temp", fasta_out_fp)
        os.remove("{0}.tmp".format(fasta_out_fp))
    os.chdir(start_dir)

def ftp_dl(x):

    # raw_path = path
    ftp_path = "ftp.ncbi.nlm.nih.gov"
    connection = FTP(ftp_path, timeout=10000)
    connection.login()

    while x:
        try:
            fp_path = x.pop()
        except IndexError:
            break
        raw_path = fp_path[0]
        fp_path = fp_path[1]
        try:
            outpath = os.path.join(raw_path, os.path.basename(fp_path))
            if os.path.isfile(outpath.strip(".gz")):
                continue

            file_size = connection.size(fp_path)
            if not os.path.isfile(outpath) or file_size != os.path.getsize(outpath):
                with open(outpath, "wb") as out_file:
                    connection.retrbinary("RETR {0}".format(fp_path), out_file.write)
            else:
                connection.sendcmd('NOOP')
        except all_errors as e:
            connection.close()
            # print(e)
            connection = FTP(ftp_path, timeout=10000)
            connection.login()

        except:
            with lock:
                with open(os.path.join(os.path.dirname(os.path.dirname(fp_path)),
                                       "/artifacts/ftp_dl.log" ),"a" ) as out:
                    out.write("{0}\n".format(fp_path))

    try:
        connection.quit()
    except:
        pass

def pull(path="",thread_count=1,databases ={"genbank"} ):
    if not path:
        string_date = datetime.datetime.now().strftime("%b-%d-%Y")
    else:
        string_date = path
    raw_path = "{}/".format(string_date)
    config_path = "artifacts/exclude.json"
    ftp_path = "ftp.ncbi.nlm.nih.gov"
    genbank_dir ="/genbank/"
    assembly_gb = "/genomes/genbank/"
    assembly_rs = "/genomes/refseq/"
    assembly_gb_summary = "assembly_summary_genbank.txt"
    assembly_rs_summary = "assembly_summary_refseq.txt"
    exclude = "suffix_exclude"

    try:
        os.makedirs(os.path.join(raw_path,"artifacts/"))
    except:
        pass
    try:
        os.makedirs(os.path.join(raw_path,"flat_files/"))
    except:
        pass

    try:
        configurations = parse_json(os.path.join(raw_path, config_path))
    except:
        configurations = {}
        configurations[exclude] = ["gbenv", "gbsyn",
                                   "gbchg" , "gbcon" , "gbnew",
                                   "gbrel", "gbtsa"]
        with open(os.path.join(raw_path,config_path), "w") as file:
            json.dump(configurations, file, sort_keys=True)
        configurations = parse_json(os.path.join(raw_path, config_path))

    exclude = set(configurations[exclude])

    connection = FTP(ftp_path)
    connection.login()
    gb_download = []
    to_download = []
    level2path = {}

    if "genbank" in databases:
        for fp in connection.nlst(genbank_dir):
            base_fp = os.path.basename(fp)
            for ind,char in enumerate(base_fp):
                try:
                    # print(ind)
                    int(char)
                    if base_fp[:ind] in exclude:
                        break
                    else:
                        gb_download.append(fp)
                        to_download.append((os.path.join(raw_path,"flat_files/"),fp))
                        break
                except:
                    continue
    level2path[b'genbank'] = gb_download


    reader = BytesIO()
    connection.retrbinary("RETR {0}{1}".format(assembly_rs,assembly_rs_summary) ,reader.write)
    reader.seek(0)

    for line in reader:
        if chr(line[0]) == "#":
            continue
        line = line.strip().split(b'\t')
        # try:
        if line[13] == b"Partial":
            continue
        try:
            if line[20].strip():
                continue
        except:
            pass
        db = line[11].strip().decode().lower().replace(" ", "_")
        if db in databases:
            try:
                temp = line[19].split(ftp_path.encode(),1)[1].decode()
                temp_path = "{0}/{1}_genomic.gbff.gz".format(temp, os.path.basename(temp))
            except:
                continue
            try:
                level2path[db.encode()].append(temp_path)
            except:
                level2path[db.encode()] = [temp_path]
            to_download.append((os.path.join(raw_path,"flat_files"),temp_path))

    reader = BytesIO()
    connection.retrbinary("RETR {0}{1}".format(assembly_gb,assembly_gb_summary) ,reader.write)
    reader.seek(0)

    for line in reader:
        if chr(line[0]) == "#":
            continue
        line = line.strip().split(b'\t')
        if line[13] == b"Partial":
            continue
        try:
            if line[20].strip():
                continue
        except:
            pass
        db = line[11].strip().decode().lower().replace(" ", "_")
        if db in databases:
            try:
                temp = line[19].split(ftp_path.encode(),1)[1].decode()
                temp_path = "{0}/{1}_genomic.gbff.gz".format(temp, os.path.basename(temp))
            except:
                continue
            try:
                level2path[db.encode()].append(temp_path)
            except:
                level2path[db.encode()] = [temp_path]

            to_download.append((os.path.join(raw_path,"flat_files/"), temp_path) )
    artifacts = [(os.path.join(raw_path,"artifacts/"),"/pub/taxonomy/taxdump.tar.gz")]
    tax_path = "/pub/taxonomy/accession2taxid/"
    for file in connection.nlst(tax_path):
        if not fnmatch.fnmatch(os.path.basename(file), 'dead*') and not fnmatch.fnmatch(file, '*md5'):
            artifacts.append((os.path.join(raw_path, "artifacts"), file))

    man = Manager()
    connection.quit()
    to_download += artifacts
    to_download = man.list(to_download)

    proc = [Process(target = ftp_dl, args=(to_download,)) for i in range(thread_count)]
    with open(os.path.join(raw_path,"artifacts", "ftp_dl.log"), "w"):
        pass
    for p in proc:
        p.start()
    for p in proc:
        p.join()

    for i in level2path.keys():
        if i.decode() in exclude:
            continue
        fp = "{0}_ff.txt".format(i.decode().replace(" ","_"))
        with open(os.path.join(raw_path,"artifacts/",fp), "w") as out_file:
            for line in level2path[i]:
                out_file.write("{0}\n".format(os.path.join(os.path.abspath(raw_path),"flat_files/", os.path.basename(line))))
    return string_date

if __name__ =="__main__":
    parser = argparse.ArgumentParser(description="TaxClipper is intended to be used to parse sequences based on NCBI taxid")

    group = parser.add_mutually_exclusive_group(required=True)
    # group.add_argument("-oc", "--oneclick", "-oneclick",action='store_true')
    # group.add_argument("-ocdl", "--oneclickdl", "-oneclickdl",action='store_true')
    # group.add_argument("-ocdc", "--oneclickdc", "-oneclickdc",action='store_true')

    group.add_argument("-pl", "--pull", "-pull",action='store_true')
    group.add_argument("-gc", "--generate-config","-generate-config", action='store_true',
                       help="generates a configuration file in current directory or as specified by output")
    group.add_argument("-bdb", "--build-database", "-build-database", action='store_true',
                       help="Builds a sequence database from a file list of NCBI flatfiles")
    # group.add_argument("-bigi", "--build-index-gi","-build-index-gi",action='store_true',
    #                    help="Builds a serialization of the fasta file using gi2taxid information")
    group.add_argument("-biacc", "--build-index-acc","-build-index-acc",action='store_true',
                       help="Builds a serialization of the fasta file using acc2taxid information")

    group.add_argument("-c", "--clip", "-clip",action='store_true')
    # group.add_argument("-ce", "--clip-exclusive","-clip-exclusive",action='store_true')
    # group.add_argument("-ci", "--clip-inclusive","-clip-inclusive",action='store_true')


    group = parser.add_argument_group()
    group.add_argument("-u", "--update", "-update", action='store_true')

    group.add_argument("-cp", "--configuration-path", "-configuration-path",
                       help="Path to configuration JSON built using generate-config flag")
    group.add_argument("-tp","--taxdump-path", "-taxdump-path", help="Path to taxdump.tar.gz file"
                       )
    group.add_argument("-sp","--serialization-path", "-serialization-path",
                       help="Path to *.p file built using build-index_* flag")
    # group.add_argument("-kp", "--keyword-path","-keyword-path",
    #                    help="Path to keyword serialization")
    group.add_argument("-fp","--fasta-path", "-fasta-path" ,
                       help="Path to fasta DB built using build-database flag")
    # group.add_argument("-g2t", "--gi-to-taxid-paths", "-gi-to-taxid-paths", nargs='*',
    #                    help="Path to gi to nucl files")
    group.add_argument("-fl","--file-list", "-file-list" ,
                       help="Path to file list of paths to GenBank Flat Files")


    group.add_argument("-a2t", "--acc-to-taxid-paths", "-acc-to-taxid-paths", nargs='*',
                       help="Path to accesion to nucl file")

    group.add_argument("-txi","--tax-id-include","-tax-id-include", nargs='*',
                       help="NCBI TaxIDs or Path to file list of TaxIDs to include")
    group.add_argument("-txe","--tax-id-exclude","-tax-id-exclude", nargs='*',
                       help="NCBI TaxIDs or Path to file list of TaxIDs to include")
    group.add_argument("-min","--minimum-length","-minimum-length", type=int,
                       help="Integer for minimum length of sequences to include")
    group.add_argument("-max","--maximum-length","-maximum-length", type=int,
                       help="Integer for maximum length of sequences to include")
    group.add_argument("-rur","--rollup-rank","-rollup-rank",
                       help="NCBI rank to set sequence headers i.e. species, genus, family et cetera")

    group.add_argument("-t", "--threads", "-threads", type=int,
                       help="Specify total threads to spawn in DB creation")

    # group.add_argument("-ow", "--overwrite", "-overwrite", help="Specify total threads to spawn in FM-index creation",
    #                    action='store_true')

    group.add_argument("-o", "--output","-output",
                       help="path for output file without extension relevant extension will be appended")

    group.add_argument("-p","--path","-path", help="Path to dated folder containing artifacts")

    args = parser.parse_args()
    if args.oneclick:
        databases = {"genbank", "complete_genome", "scaffold", "contig", "chromosome"}
        exclude = {"complete_Genome", "contig"}
        if args.threads:
            threads = args.threads
        else:
            threads = 1
        if args.path:
            dl_folder = args.path
        else:
            dl_folder = ""
        dl_folder = pull(thread_count=threads,excluded=exclude, path=dl_folder)

        for fp in os.listdir(os.path.join(dl_folder,"artifacts/")):
            if fnmatch.fnmatch(fp, "*_ff.txt"):
                build_db(os.path.join(dl_folder,"artifacts/",fp), os.path.join(dl_folder,"artifacts/","{0}.fas".format(fp.split("_ff")[0])),os.devnull, os.devnull, threads, os.devnull)

        arguments = oneclickjson(dl_folder)
        with Pool(threads) as p:
            p.starmap(acc_serialization, [(argument['acc-to-taxid-paths'], argument['fasta-path'], argument['taxdump-path']) for argument in arguments ] )
                      # (arguments['acc-to-taxid-paths'], arguments['fasta-path'], arguments['taxdump-path'])
#        shutil.rmtree(os.path.join(dl_folder,"flat_files"))

    elif args.pull:
        if args.threads and args.tax_id_exclude:
            pull(thread_count=args.threads, excluded=set(args.tax_id_exclude))
        elif args.threads:
            pull(thread_count=args.threads)
        elif args.tax_id_exclude:
            pull(excluded=set(args.tax_id_exclude))
        else:
            pull()
    else:

        if args.configuration_path:
            arguments = arg_unwrappers(args,parse_json(args))
        else:
            arguments = arg_unwrappers(args)

        if args.build_database:
            if args.threads:
                threads = args.threads
            else:
                threads = 1
            if args.file_list:
            # build_db("chromosomes.list","chromosome.fasta", "test.kw", "test.src", 16, "test.g2w")
            # build_db("complete_genome.list","complete.fasta", "test.kw", "test.src", 16, "test.g2w")
                build_db(*"{2} {0} {0}.kw {0}.src {1} {0}.g2w".format(args.output, threads, args.file_list).split())
            else:
                print("FASTA Database creation requires a path to a file list of GenBank Flat Files")
        elif args.build_index_acc:
            if arguments['taxdump-path'] and arguments['fasta-path'] and arguments['acc-to-taxid-paths']:
                acc_serialization(arguments['acc-to-taxid-paths'],arguments['fasta-path'],arguments['taxdump-path'])
            else:
                parser.error("Serialization requires paths to taxdump, fasta database and accession2taxid files")
            if args.update:
                arguments['serialization-path'] = os.path.abspath(args.output+".p")
                gen_json(arguments, args)

        elif args.clip:
            clip(args.tax_id_include,arguments['rollup-rank'], args.tax_id_exclude,
                 args.output,arguments['minimum-length'], arguments['maximum-length'], arguments['fasta-path'], arguments['serialization-path'])
        elif args.update and args.configuration_path:
            gen_json(arg_unwrappers(args,parse_json(args)),args)

        elif args.generate_config:
            gen_json(arg_unwrappers(args),args)
