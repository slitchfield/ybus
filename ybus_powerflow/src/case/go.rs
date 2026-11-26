//! Reading the input files of ARPA-E's Grid Optimization (GO) Competition Challenge 1, the
//! security-constrained AC optimal power flow challenge (`docs/GoParserChecklist.md`, goal 4).
//!
//! An instance is four ASCII files: `case.raw` (the network, in a PSS/E revision 33 subset),
//! `case.rop` (generator costs), `case.inl` (participation factors) and `case.con`
//! (contingencies). The Challenge 1 problem formulation's Appendix A defines them, and its
//! Appendix C how they're read. Errors name the file and the line, as
//! [`PowerflowError::GoParse`](crate::PowerflowError::GoParse).

use std::fmt;

mod text;

/// One of the four files of a GO Challenge 1 instance, as errors name it.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
#[non_exhaustive]
pub enum GoFile {
    /// `case.raw`: the network, in PSS/E's raw format (revision 33).
    Raw,
    /// `case.rop`: the generator cost curves, in PSS/E's ROP format.
    Rop,
    /// `case.inl`: the participation factors, in PSS/E's INLF format.
    Inl,
    /// `case.con`: the contingencies, in PSS/E's contingency description format.
    Con,
}

impl GoFile {
    /// The file's name in an instance's folders, such as `case.raw`.
    pub fn name(self) -> &'static str {
        match self {
            GoFile::Raw => "case.raw",
            GoFile::Rop => "case.rop",
            GoFile::Inl => "case.inl",
            GoFile::Con => "case.con",
        }
    }
}

impl fmt::Display for GoFile {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.name())
    }
}
