use noodles::bam;
use noodles::sam::alignment::Record;
use noodles::sam::record::ReadName;
use crate::bed::{BED, Score, Strand};
use indexmap::map::IndexMap;
use noodles::sam::header::ReferenceSequence;
use std::collections::hash_map::RandomState;
use serde::{Serialize, Serializer, Deserialize, Deserializer};
use extsort::sorter::Sortable;

pub struct BamRecord(Record);

/// An alignment record.
#[derive(Clone, Debug, PartialEq)]
pub struct SamRecord {
    read_name: Option<String>,
    flags: u16,
    reference_sequence_id: Option<usize>,
    alignment_start: Option<Position>,
    mapping_quality: Option<MappingQuality>,
    cigar: Cigar,
    mate_reference_sequence_id: Option<usize>,
    mate_alignment_start: Option<Position>,
    template_length: i32,
    sequence: Sequence,
    quality_scores: QualityScores,
    data: Data,
}

impl SamRecord {
    pub fn from_record(rec: Record) -> Self {
        SamRecord {
            read_name: rec.read_name().map(|x| x.to_string()),
            flags: rec.flags().bits()

        }

    }

    pub fn to_record()
}

impl Serialize for BamRecord {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        self.0.read_name().map(|x| <ReadName as AsRef<[u8]>>::as_ref(x)).serialize(serializer)?;
        todo!()
    }
}

impl<'de> Deserialize<'de> for BamRecord {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        todo!()
    }
}

impl Sortable for BamRecord {
    fn encode<W: std::io::Write>(&self, writer: &mut W) {
        bincode::serialize_into(writer, self).unwrap();
    }

    fn decode<R: std::io::Read>(reader: &mut R) -> Option<Self> {
        bincode::deserialize_from(reader).ok()
    }
}


/*
pub fn bam_to_bed(refs: &IndexMap<String, ReferenceSequence, RandomState>,
                  record: &Record)
                  -> Option<BED<6>> {
    let chr = record.reference_sequence(refs)?.unwrap().name().as_str().to_string();
    let start_loc: i32 = record.alignment_start()?.into();
    let end_loc: i32 = record.alignment_end()?.unwrap().into();
    let name = record.read_name().unwrap().to_str().unwrap().to_string();
    let score = record.mapping_quality().map(|x| {
        let x_: u8 = x.into();
        Score::try_from(x_ as u16).unwrap()
    });
    let strand = if record.flags().is_reverse_complemented()
        { Strand::Reverse } else { Strand::Forward };
    Some(BED::new(chr, (start_loc - 1) as u64, end_loc as u64,
                  Some(name), score, Some(strand), Default::default()))
}
*/