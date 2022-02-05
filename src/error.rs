use thiserror::Error;

pub type Result<T> = std::result::Result<T, Error>;

#[derive(Debug, Error)]
pub enum Error {
    #[error("Decimal conversion syntax error")]
    ConversionSyntax,
    #[error("Decimal division by zero")]
    DivisionByZero,
    #[error("Decimal overflow")]
    Overflow,
    #[error("Exceeds range of conversion target")]
    ExceedsConversionTargetRange,
}
