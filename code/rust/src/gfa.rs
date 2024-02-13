use std::boxed::Box;
use std::collections::HashMap;
use std::path::Path;
use std::str::FromStr;

pub enum Category {
    Segment,
}

impl std::fmt::Display for Category {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}",
            match self {
                Self::Segment => "S",
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
///```
pub struct Line {
    pub category: Category,
    pub name: String,
    pub sequence: String,
    pub tags: Vec<String>,
}

impl std::str::FromStr for Line {
    type Err = Error;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let mut toks = s.split_terminator('\t');
        let category =
            Category::from_str(toks.next().ok_or_else(|| {
                Error(format!("Expected token type in malformatted string {s}"))
            })?)?;
        let name = toks
            .next()
            .ok_or_else(|| Error(format!("Expected name type in malformatted string {s}")))?
            .to_owned();
        let sequence = toks
            .next()
            .ok_or_else(|| Error(format!("Expected seq type in malformatted string {s}")))?
            .to_owned();
        let tags = toks.map(|x| x.to_string()).collect::<Vec<_>>();
        Ok(Self {
            category,
            name,
            sequence,
            tags,
        })
    }
}

impl std::fmt::Display for Line {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let res = write!(f, "{}\t{}\t{}", self.category, self.name, self.sequence);
        res.and_then(|x| {
            if !self.tags.is_empty() {
                write!(
                    f,
                    "\t{}",
                    itertools::Itertools::intersperse(self.tags.iter().map(|x| &x[..]), "\t")
                        .collect::<String>()
                )
            } else {
                Ok(x)
            }
        })
    }
}

pub struct File {
    pub lines: Vec<Line>,
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
        let mut lines = Vec::with_capacity(8192);
        let mut name2idx = HashMap::new();
        for (idx, input) in reader.lines().enumerate() {
            let input = input?;
            let line = Line::from_str(&input)?;
            name2idx.insert(line.name.clone(), idx);
            lines.push(line);
        }
        Ok(Self { lines, name2idx })
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
}
