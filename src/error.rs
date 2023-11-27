#[derive(PartialEq, Debug)]
pub enum Error
{
    SampleTooSmall,
    IndexOutOfRange,
    WrongInputSize,
    ZeroDenominator,
}
