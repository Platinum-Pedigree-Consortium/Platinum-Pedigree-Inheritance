use std::boxed::Box;
use std::collections::{BTreeMap, HashMap};
use std::path::Path;
use std::str::FromStr;

#[repr(u8)]
#[derive(Clone, Copy)]
pub enum Category {
    Segment,
    Link,
}

impl std::default::Default for Category {
    fn default() -> Self {
        Self::Segment
    }
}

impl std::fmt::Display for Category {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}",
            match self {
                Self::Segment => "S",
                Self::Link => "L",
            }
        )
    }
}

#[derive(Debug, PartialEq, Eq)]
pub struct Error(pub String);

impl std::fmt::Display for Error {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "GFAError{{msg: {}}}", self.0)
    }
}

impl std::error::Error for Error {}

impl std::str::FromStr for Category {
    type Err = Error;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "S" => Ok(Category::Segment),
            "L" => Ok(Category::Link),
            _ => Err(Error(format!("Failed to get category: {s}"))),
        }
    }
}

///
///
///```
/// use concordance::gfa::Line;
/// use std::str::FromStr;
/// let line = Line::from_str("S\tS1\tACGT\tAC:Z:0").expect("Failed to parse line");
/// eprintln!("Line: {line}");
/// assert_eq!(line.to_string(), "S\tS1\tACGT\tAC:Z:0");
/// assert_eq!(line.name(), "S1");
///```
#[derive(Clone, Default)]
pub struct Line {
    pub category: Category,
    pub name: Option<String>,
    pub sequence: Option<String>,
    pub tags: Vec<String>,
    pub from: Option<String>,
    pub to: Option<String>,
    pub cigar: Option<String>,
    pub from_orientation: bool,
    pub to_orientation: bool,
}

impl Line {
    #[must_use]
    pub fn name(&self) -> String {
        if let Some(name) = &self.name {
            name.clone()
        } else {
            format!("{:?}->{:?}", self.from().unwrap(), self.to().unwrap())
        }
    }

    #[must_use]
    fn from(&self) -> Option<&String> {
        self.from.as_ref()
    }

    #[must_use]
    fn to(&self) -> Option<&String> {
        self.to.as_ref()
    }
}

/// Store orientation as a bool. True = forward strand.
#[must_use]
fn orientation(s: &str) -> bool {
    match s {
        "+" => true,
        "-" => false,
        _ => {
            panic!("Invalid orientation {s}");
        }
    }
}

#[must_use]
fn fmtstrand(x: bool) -> char {
    if x {
        '+'
    } else {
        '-'
    }
}

impl std::str::FromStr for Line {
    type Err = Error;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let mut toks = s.split_terminator('\t');
        let category =
            Category::from_str(toks.next().ok_or_else(|| {
                Error(format!("Expected token type in malformatted string {s}"))
            })?)?;
        match category {
            Category::Segment => {
                let name = Some(
                    toks.next()
                        .ok_or_else(|| {
                            Error(format!("Expected name type in malformatted string {s}"))
                        })?
                        .to_owned(),
                );
                let sequence = Some(
                    toks.next()
                        .ok_or_else(|| {
                            Error(format!("Expected seq type in malformatted string {s}"))
                        })?
                        .to_owned(),
                );
                let tags = toks
                    .map(std::string::ToString::to_string)
                    .collect::<Vec<_>>();
                Ok(Self {
                    category,
                    name,
                    sequence,
                    tags,
                    ..Default::default()
                })
            }
            Category::Link => {
                let from = Some(
                    toks.next()
                        .ok_or_else(|| Error(format!("Expected from segment id in s {s}")))
                        .map(|x| x.to_owned())?,
                );
                let from_orientation = orientation(
                    toks.next()
                        .ok_or_else(|| Error("Expected from orientation".to_string()))?,
                );
                let to = Some(
                    toks.next()
                        .ok_or_else(|| Error(format!("Expected to segment id in s {s}")))
                        .map(|x| x.to_owned())?,
                );
                let to_orientation = orientation(
                    toks.next()
                        .ok_or_else(|| Error("Expected to orientation".to_string()))?,
                );
                let cigar = Some(
                    toks.next()
                        .map(|x| x.to_owned())
                        .ok_or_else(|| Error("Expected cigar".into()))?,
                );
                let tags = toks
                    .map(std::string::ToString::to_string)
                    .collect::<Vec<_>>();
                Ok(Self {
                    category,
                    tags,
                    cigar,
                    from,
                    to,
                    from_orientation,
                    to_orientation,
                    ..Default::default()
                })
            }
        }
    }
}

impl std::fmt::Display for Line {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self.category {
            Category::Segment => {
                let res = write!(
                    f,
                    "{}\t{}\t{}",
                    self.category,
                    self.name.as_ref().unwrap(),
                    self.sequence.as_ref().unwrap()
                );
                res.and_then(|x| {
                    if self.tags.is_empty() {
                        Ok(x)
                    } else {
                        write!(
                            f,
                            "\t{}",
                            itertools::Itertools::intersperse(
                                self.tags.iter().map(|x| &x[..]),
                                "\t"
                            )
                            .collect::<String>()
                        )
                    }
                })
            }
            Category::Link => {
                let res = write!(
                    f,
                    "{}\t{}\t{}\t{}\t{}\t{}",
                    self.category,
                    self.from().unwrap(),
                    fmtstrand(self.from_orientation),
                    self.to().unwrap(),
                    fmtstrand(self.to_orientation),
                    self.cigar.as_ref().unwrap()
                );
                res.and_then(|x| {
                    if self.tags.is_empty() {
                        Ok(x)
                    } else {
                        write!(
                            f,
                            "\t{}",
                            itertools::Itertools::intersperse(
                                self.tags.iter().map(|x| &x[..]),
                                "\t"
                            )
                            .collect::<String>()
                        )
                    }
                })
            }
        }
    }
}

/// Holds segments + links + edges.
/// Look up the numeric id for a segment with `name2idx`.
/// Find excident edges in `edges` using a `BTreeMap` range query.
pub struct File {
    pub segments: Vec<Line>,
    pub links: Vec<Line>,
    pub edges: BTreeMap<(usize, usize), usize>,
    pub name2idx: HashMap<String, usize>,
}

impl File {
    /// Build `gfa::File` from a type implementing `std::io::Read`
    ///
    /// # Errors
    /// 1. Lines not valid UTF-8.
    /// 2. Lines are not valid gfa.
    /// 3. Lines are gfa lines we don't support yet.
    pub fn from_reader(x: impl std::io::Read) -> Result<Self, Box<dyn std::error::Error>> {
        use std::io::BufRead;
        let reader = std::io::BufReader::new(x);
        let mut segments = Vec::with_capacity(8192);
        let mut links = Vec::with_capacity(8192);
        let mut edges = BTreeMap::<(usize, usize), usize>::new();
        let mut name2idx = HashMap::new();
        let mut line_num = 0usize;
        for input in reader.lines() {
            line_num += 1;
            if line_num % 10000 == 0 {
                log::info!("Parsed {line_num} lines from gfa file");
            }
            let input = input?;
            let line = Line::from_str(&input)?;
            match line.category {
                Category::Link => {
                    let to = line.to().and_then(|x| name2idx.get(x)).expect("Missing to");
                    let from = line
                        .from()
                        .and_then(|x| name2idx.get(x))
                        .expect("Missing from");
                    edges.insert((*to, *from), links.len());
                    links.push(line);
                }
                Category::Segment => {
                    name2idx.insert(line.name.as_ref().unwrap().clone(), segments.len());
                    segments.push(line);
                }
            }
        }
        Ok(Self {
            segments,
            links,
            edges,
            name2idx,
        })
    }

    /// Build `gfa::File` from a `Path`.
    ///
    /// # Errors
    /// 1. File does not exist.
    /// 2. Cannot read from file.
    /// 3. Lines not valid UTF-8.
    /// 4. Lines are not valid gfa.
    /// 5. Lines are gfa lines we don't support yet.
    pub fn from_path(x: impl AsRef<Path>) -> Result<Self, Box<dyn std::error::Error>> {
        let path = std::fs::File::open(x.as_ref())?;
        Self::from_reader(path)
    }

    /// Returns all excident edges for a given segment id.
    /// # Panics
    /// Panics if start > end, which can never happen in software as `0 < usize::MAX`.
    #[must_use]
    pub fn out_edges(
        &self,
        idx: usize,
    ) -> std::collections::btree_map::Range<'_, (usize, usize), usize> {
        self.edges.range((
            std::ops::Bound::Included(&(idx, 0)),
            std::ops::Bound::Included(&(idx, usize::MAX)),
        ))
    }
}

impl std::fmt::Display for File {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "GFAFile{{{} segments, {} links}}",
            self.segments.len(),
            self.links.len()
        )
    }
}
