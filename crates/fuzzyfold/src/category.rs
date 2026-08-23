#[derive(Clone, Eq, Hash, PartialEq)]
pub enum Category {
    Within(usize),
    Between(usize, usize),
    WithRest(usize),
}

impl Category {
    pub fn to_key(&self) -> String {
        match self {
            Category::Within(a)     => format!("Within_{}", a + 1),
            Category::Between(a, b) => format!("Between_{}_{}", a + 1, b + 1),
            Category::WithRest(a)   => format!("WithRest_{}", a + 1),
        }
    }

    pub fn from_key(s: &str) -> anyhow::Result<Self> {
        let parts: Vec<&str> = s.split('_').collect();
        match parts.as_slice() {
            ["Within", a]     => Ok(Category::Within(a.parse::<usize>()? - 1)),
            ["Between", a, b] => Ok(Category::Between(a.parse::<usize>()? - 1, b.parse::<usize>()? - 1)),
            ["WithRest", a]   => Ok(Category::WithRest(a.parse::<usize>()? - 1)),
            _ => anyhow::bail!("unknown category key: {}", s),
        }
    }
}