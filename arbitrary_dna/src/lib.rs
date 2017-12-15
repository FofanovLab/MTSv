extern crate quickcheck;

use quickcheck::{Arbitrary, Gen};

pub type Dna5Sequence = Vec<Dna5Base>;

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct Dna5Base(pub u8);

impl Arbitrary for Dna5Base {
    fn arbitrary<G: Gen>(g: &mut G) -> Self {
        let rand = g.next_u32();

        Dna5Base(match rand % 5 {
            0 => b'A',
            1 => b'C',
            2 => b'G',
            3 => b'T',
            4 => b'N',
            _ => unreachable!(),
        })
    }

    fn shrink(&self) -> Box<Iterator<Item = Self>> {
        Dna5BaseShrinker::new(*self)
    }
}

struct Dna5BaseShrinker {
    current: Dna5Base,
}

impl Dna5BaseShrinker {
    pub fn new(x: Dna5Base) -> Box<Iterator<Item=Dna5Base>> {
        Box::new(Dna5BaseShrinker { current: x })
    }
}

impl Iterator for Dna5BaseShrinker {
    type Item = Dna5Base;
    fn next(&mut self) -> Option<Dna5Base> {
        let to_ret = match self.current.0 {
            b'N' => Some(Dna5Base(b'T')),
            b'T' => Some(Dna5Base(b'G')),
            b'G' => Some(Dna5Base(b'C')),
            b'C' => Some(Dna5Base(b'A')),
            _ => None,
        };

        if let Some(b) = to_ret {
            self.current = b;
        }

        to_ret
    }
}
