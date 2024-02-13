use std::boxed::Box;
use std::str::FromStr;

pub enum Category {
    Segment,
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

pub struct Line {
    pub category: Category,
    pub name: String,
    pub sequence: Vec<u8>,
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
            .as_bytes()
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

pub struct File {
    pub lines: Vec<Line>,
}

impl File {
    pub fn from_reader(x: impl std::io::Read) -> Result<Self, Box<dyn std::error::Error>> {
        use std::io::BufRead;
        let reader = std::io::BufReader::new(x);
        let mut lines = Vec::with_capacity(8192);
        for input in reader.lines() {
            let input = input?;
            lines.push(Line::from_str(&input)?);
        }
        Ok(Self { lines })
    }
}
