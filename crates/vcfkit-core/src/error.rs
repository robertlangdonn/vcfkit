use thiserror::Error;

#[derive(Debug, Error)]
pub enum VcfkitError {
    #[error("I/O error: {0}")]
    Io(#[from] std::io::Error),
    #[error("{0}")]
    Other(String),
}
