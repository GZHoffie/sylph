//Various byte-tables and hashing methods are taken from miniprot by Heng Li. Attached below is their license:
//The MIT License

// **** miniprot LICENSE ***
//Copyright (c) 2022-     Dana-Farber Cancer Institute
//
//Permission is hereby granted, free of charge, to any person obtaining
//a copy of this software and associated documentation files (the
//"Software"), to deal in the Software without restriction, including
//without limitation the rights to use, copy, modify, merge, publish,
//distribute, sublicense, and/or sell copies of the Software, and to
//permit persons to whom the Software is furnished to do so, subject to
//the following conditions:
//
//The above copyright notice and this permission notice shall be
//included in all copies or substantial portions of the Software.
//
//THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
//EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
//MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
//NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
//BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
//ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
//CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
//SOFTWARE.
//******************************


use std::collections::HashMap;

// bytecheck can be used to validate your data if you want
use std::hash::{BuildHasherDefault, Hasher};
use std::collections::HashSet;
use smallvec::SmallVec;
use serde::{Deserialize, Serialize, Serializer, Deserializer, de::Visitor};
use fxhash::FxHashMap;

#[derive(Deserialize, Serialize, Debug, PartialEq)]
pub enum AdjustStatus {
    Lambda(f64),
    Low,
    High,
}

impl Default for AdjustStatus {
    fn default() -> Self {AdjustStatus::Low }
}

pub type Kmer = u64;
pub const BYTE_TO_SEQ: [u8; 256] = [
    0, 1, 2, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
];

pub const BYTE_TO_ERROR_RATE: [f64; 256] = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.7943282347242815, 0.6309573444801932, 0.5011872336272722, 0.3981071705534972, 0.31622776601683794, 0.251188643150958, 0.19952623149688797, 0.15848931924611134, 0.12589254117941673, 0.1, 0.07943282347242814, 0.06309573444801933, 0.05011872336272722, 0.039810717055349734, 0.03162277660168379, 0.025118864315095794, 0.0199526231496888, 0.015848931924611134, 0.012589254117941675, 0.01, 0.007943282347242814, 0.00630957344480193, 0.005011872336272725, 0.003981071705534973, 0.0031622776601683794, 0.0025118864315095794, 0.001995262314968879, 0.001584893192461114, 0.0012589254117941675, 0.001, 0.0007943282347242813, 0.000630957344480193, 0.0005011872336272725, 0.00039810717055349735, 0.00031622776601683794, 0.00025118864315095795, 0.00019952623149688788, 0.00015848931924611142, 0.00012589254117941674, 0.0001, 7.943282347242822e-05, 6.309573444801929e-05, 5.011872336272725e-05, 3.9810717055349695e-05, 3.1622776601683795e-05, 2.5118864315095822e-05, 1.9952623149688786e-05, 1.584893192461114e-05, 1.2589254117941661e-05, 1e-05, 7.943282347242822e-06, 6.30957344480193e-06, 5.011872336272725e-06, 3.981071705534969e-06, 3.162277660168379e-06, 2.5118864315095823e-06, 1.9952623149688787e-06, 1.584893192461114e-06, 1.2589254117941661e-06, 1e-06, 7.943282347242822e-07, 6.30957344480193e-07, 5.011872336272725e-07, 3.981071705534969e-07, 3.162277660168379e-07, 2.5118864315095823e-07, 1.9952623149688787e-07, 1.584893192461114e-07, 1.2589254117941662e-07, 1e-07, 7.943282347242822e-08, 6.30957344480193e-08, 5.011872336272725e-08, 3.981071705534969e-08, 3.162277660168379e-08, 2.511886431509582e-08, 1.9952623149688786e-08, 1.5848931924611143e-08, 1.2589254117941661e-08, 1e-08, 7.943282347242822e-09, 6.309573444801943e-09, 5.011872336272715e-09, 3.981071705534969e-09, 3.1622776601683795e-09, 2.511886431509582e-09, 1.9952623149688828e-09, 1.584893192461111e-09, 1.2589254117941663e-09, 1e-09, 7.943282347242822e-10, 6.309573444801942e-10, 5.011872336272714e-10, 3.9810717055349694e-10, 3.1622776601683795e-10, 2.511886431509582e-10, 1.9952623149688828e-10, 1.584893192461111e-10, 1.2589254117941662e-10, 1e-10, 7.943282347242822e-11, 6.309573444801942e-11, 5.011872336272715e-11, 3.9810717055349695e-11, 3.1622776601683794e-11, 2.5118864315095823e-11, 1.9952623149688828e-11, 1.5848931924611107e-11, 1.2589254117941662e-11, 1e-11, 7.943282347242821e-12, 6.309573444801943e-12, 5.011872336272715e-12, 3.9810717055349695e-12, 3.1622776601683794e-12, 2.5118864315095823e-12, 1.9952623149688827e-12, 1.584893192461111e-12, 1.258925411794166e-12, 1e-12, 7.943282347242822e-13, 6.309573444801942e-13, 5.011872336272715e-13, 3.981071705534969e-13, 3.162277660168379e-13, 2.511886431509582e-13, 1.9952623149688827e-13, 1.584893192461111e-13, 1.2589254117941663e-13, 1e-13, 7.943282347242822e-14, 6.309573444801943e-14, 5.0118723362727144e-14, 3.9810717055349693e-14, 3.1622776601683796e-14, 2.5118864315095823e-14, 1.9952623149688828e-14, 1.584893192461111e-14, 1.2589254117941662e-14, 1e-14, 7.943282347242822e-15, 6.309573444801943e-15, 5.0118723362727146e-15, 3.9810717055349695e-15, 3.1622776601683794e-15, 2.511886431509582e-15, 1.995262314968883e-15, 1.584893192461111e-15, 1.2589254117941663e-15, 1e-15, 7.943282347242821e-16, 6.309573444801943e-16, 5.011872336272715e-16, 3.9810717055349695e-16, 3.1622776601683793e-16, 2.511886431509582e-16, 1.995262314968883e-16, 1.5848931924611109e-16, 1.2589254117941662e-16, 1e-16, 7.943282347242789e-17, 6.309573444801943e-17, 5.0118723362727144e-17, 3.9810717055349855e-17, 3.1622776601683796e-17, 2.5118864315095718e-17, 1.9952623149688827e-17, 1.584893192461111e-17, 1.2589254117941713e-17, 1e-17, 7.94328234724279e-18, 6.309573444801943e-18, 5.011872336272715e-18, 3.981071705534985e-18, 3.1622776601683795e-18, 2.5118864315095718e-18, 1.995262314968883e-18, 1.5848931924611109e-18, 1.2589254117941713e-18, 1e-18, 7.943282347242789e-19, 6.309573444801943e-19, 5.011872336272715e-19, 3.9810717055349853e-19, 3.162277660168379e-19, 2.5118864315095717e-19, 1.995262314968883e-19, 1.584893192461111e-19, 1.2589254117941713e-19, 1e-19, 7.94328234724279e-20, 6.309573444801943e-20, 5.011872336272715e-20, 3.9810717055349855e-20, 3.162277660168379e-20, 2.511886431509572e-20, 1.9952623149688828e-20, 1.5848931924611108e-20, 1.2589254117941713e-20, 1e-20, 7.943282347242789e-21, 6.309573444801943e-21, 5.011872336272714e-21, 3.981071705534986e-21, 3.1622776601683792e-21, 2.511886431509572e-21, 1.9952623149688827e-21, 1.5848931924611108e-21, 1.2589254117941713e-21, 1e-21, 7.943282347242789e-22, 6.309573444801943e-22, 5.011872336272715e-22, 3.9810717055349856e-22, 3.1622776601683793e-22, 2.511886431509572e-22, 1.9952623149688828e-22, 1.584893192461111e-22, 1.2589254117941713e-22, 1e-22, 7.943282347242789e-23, 6.309573444801943e-23];

#[inline]
pub fn mm_hash(bytes: &[u8]) -> usize {
    let mut key = usize::from_ne_bytes(bytes.try_into().unwrap()) as usize;
    key = (!key).wrapping_add(key << 21); // key = (key << 21) - key - 1;
    key = key ^ key >> 24;
    key = (key.wrapping_add(key << 3)).wrapping_add(key << 8); // key * 265
    key = key ^ key >> 14;
    key = (key.wrapping_add(key << 2)).wrapping_add(key << 4); // key * 21
    key = key ^ key >> 28;
    key = key.wrapping_add(key << 31);
    return key;
}

pub struct MMHasher {
    hash: usize,
}

impl Hasher for MMHasher {
    #[inline]
    fn write(&mut self, bytes: &[u8]) {
        self.hash = mm_hash(bytes);
    }
    #[inline]
    fn finish(&self) -> u64 {
        self.hash as u64
    }
}

impl Default for MMHasher {
    #[inline]
    fn default() -> MMHasher {
        MMHasher { hash: 0 }
    }
}

//Implement minimap2 hashing, will test later.
pub type MMBuildHasher = BuildHasherDefault<MMHasher>;
pub type MMHashMap<K, V> = HashMap<K, V, MMBuildHasher>;
pub type MMHashSet<K> = HashSet<K, MMBuildHasher>;

/// `serde` helpers to improve serialization of the `FxHashMap` storing k-mer counts.
/// 
/// Encoding the `FxHashMap` as a sequence instead of a map speeds up serialize
/// and deserialize by a magnitude.
mod kmer_counts {
    use super::*;

    struct FxHashMapVisitor;
    
    impl<'a> Visitor<'a> for FxHashMapVisitor {
        type Value = FxHashMap<Kmer, u32>;

        fn expecting(&self, formatter: &mut std::fmt::Formatter) -> std::fmt::Result {
            formatter.write_str("a sequence of kmer counts")
        }

        fn visit_seq<A>(self, mut seq: A) -> Result<Self::Value, A::Error>
            where
                A: serde::de::SeqAccess<'a>
        {
            let mut counts = match seq.size_hint() {
                Some(size) => FxHashMap::with_capacity_and_hasher(size, Default::default()),
                None => FxHashMap::default(),
            };
            while let Some(item) = seq.next_element::<(Kmer, u32)>()? {
                counts.insert(item.0, item.1);
            }
            Ok(counts)
        }
    }

    pub fn serialize<S>(
        kmer_counts: &FxHashMap<Kmer, u32>, 
        serializer: S
    ) -> Result<S::Ok, S::Error> 
    where S: Serializer {
        serializer.collect_seq(kmer_counts.into_iter())
    }

    pub fn deserialize<'de, D>(deserializer: D) -> Result<FxHashMap<Kmer, u32>, D::Error> where D: Deserializer<'de> {
        deserializer.deserialize_seq(FxHashMapVisitor)
    }
}

#[derive(Default, Deserialize, Serialize, Debug, PartialEq)]
pub struct SequencesSketch{
    #[serde(with = "kmer_counts")]
    pub kmer_counts: FxHashMap<Kmer, u32>,
    pub c: usize,
    pub k: usize,
    pub file_name: String,
    pub sample_name: Option<String>,
    pub paired: bool,
    pub mean_read_length: f64,
    pub average_kmers_per_read: f64,
    pub average_error_rate: f64,
}

impl SequencesSketch{
    pub fn new(file_name: String, c: usize, k: usize, paired: bool, sample_name: Option<String>, mean_read_length: f64) -> SequencesSketch{
        return SequencesSketch{kmer_counts : HashMap::default(), file_name, c, k, paired, sample_name, mean_read_length}
    }
}

#[derive(Deserialize, Serialize, Debug, PartialEq, Hash, PartialOrd, Eq, Ord, Default, Clone)]
pub struct GenomeSketch{
    pub genome_kmers: Vec<Kmer>,
    pub pseudotax_tracked_nonused_kmers: Option<Vec<Kmer>>,
    pub file_name: String,
    pub first_contig_name: String,
    pub c: usize,
    pub k: usize,
    pub gn_size: usize,
    pub min_spacing: usize,
}

#[derive(Deserialize, Serialize, Debug, PartialEq)]
#[derive(Default, Clone)]
pub struct MultGenomeSketch{
    pub genome_kmer_index: Vec<(Kmer,SmallVec<[u32;1]>)>,
    pub file_names: Vec<String>,
    pub contig_names: Vec<String>,
    pub c: usize,
    pub k: usize,
}

#[derive(Debug, PartialEq)]
pub struct AniResult<'a>{
    pub naive_ani: f64,
    pub final_est_ani: f64,
    pub naive_ani1: f64,
    pub final_est_ani1: f64,
    pub final_est_cov: f64,
    pub seq_name: String,
    pub gn_name: &'a str,
    pub contig_name: &'a str,
    pub mean_cov: f64,
    pub median_cov: f64,
    pub containment_index: (usize,usize),
    pub conditional_containment_index: (usize,usize),
    pub lambda: AdjustStatus,
    pub ani_ci: (Option<f64>,Option<f64>),
    pub lambda_ci: (Option<f64>,Option<f64>),
    pub genome_sketch: &'a GenomeSketch,
    pub rel_abund: Option<f64>,
    pub seq_abund: Option<f64>,
    pub kmers_lost: Option<usize>,

}
