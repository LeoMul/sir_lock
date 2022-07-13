use core::panic;
use std::num::NonZeroU64;

use{
    structopt::StructOpt,
    std::{
        fs::File,
        io::{BufWriter, Write, BufReader},
    },
    crate::sir_model::*,
    super::*,
};
use net_ensembles::{sampling::{HistogramFast}};


use {
    rand_pcg::Pcg64,
    net_ensembles::sampling::Rewl,
};

#[derive(StructOpt, Debug, Clone)]
/// Replica exchange entropic sampling for high degree
pub struct ReesOpts
{
    /// Specify the save file of REWL
    #[structopt(long)]
    pub save_file: String,

    #[structopt(long)]
    /// Add 'hours' to the available time
    pub hours: Option<NonZeroU64>,

    #[structopt(long)]
    /// Add 'minutes' to the available time
    pub minutes: Option<NonZeroU64>,

    #[structopt(long)]
    /// Add 'days' to the available time
    pub days: Option<NonZeroU64>,

    /// Do not save the entropic sampling state at the end
    #[structopt(long)]
    pub no_save: bool,

    #[structopt(long)]
    /// How many curves should be printed (per walker)
    pub print_count: NonZeroU64

}

impl ReesOpts
{
    pub fn allowed_seconds(&self) -> u64
    {
        let mut secs = 0;
        if let Some(min) = self.minutes
        {
            secs += min.get() * 60;
        }

        if let Some(hours) = self.hours
        {
            secs += hours.get() * 60 * 60;
        }

        if let Some(days) = self.days
        {
            secs += days.get() * 60 * 60 * 24;
        }

        if secs == 0 {
            panic!("Specified time is 0. Unable to process anything without time")
        }

        secs
    }
    pub fn execute(self, instant: std::time::Instant)
    {
        execute_entropic_sampling(self,instant)
    }
}

pub struct LogfilePrinter
{
    printer: BufWriter<File>
}

impl LogfilePrinter
{
    pub fn new(name: &str) -> Self
    {
        let file = File::create(name)
            .unwrap();
        let buf = BufWriter::new(file);
        
        Self{
            printer: buf
        }
    }
}

impl Write for LogfilePrinter
{
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize>
    {
        std::io::stdout().write_all(buf)?;
        self.printer.write(buf)
    }

    fn flush(&mut self) -> std::io::Result<()>
    {
        self.printer.flush()
    }
}


pub type BALDRewl = Rewl<BALargeDeviationWithLocks, Pcg64, HistogramFast<u32>, u32, MarkovStepWithLocks, ()>;


pub fn deserialize_from_file(filename: &str) -> (BALDRewl, Vec<String>)
{
    let file = File::open(filename).expect("unable to open save file");
    let reader = BufReader::new(file);

    let res: Result<
        (
            BALDRewl, 
            String
        ), _> =  bincode::deserialize_from(reader);

        match res
        {
            Ok((rewl, json_string)) => {
                (
                    rewl,
                    vec![json_string]
                )
            },
            _ => 
            {
                let file = File::open(filename).expect("unable to open save file");
                let reader = BufReader::new(file);
                bincode::deserialize_from(reader).expect("unable to parse bincode file")
            }
        }
}