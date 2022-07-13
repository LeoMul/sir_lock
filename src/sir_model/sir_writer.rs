use std::fmt::Display;

use serde_json::Value;
use std::mem::ManuallyDrop;
use std::ops::DerefMut;
use std::process::Command;

use{
    std::{
        fs::File,
        io::{Write, BufWriter}
    },
    net_ensembles::{GenericGraph, AdjContainer},
    super::*
};

pub type CurveWriter = BufWriter<File>;
pub struct SirWriter
{
    pub writer_s: ManuallyDrop<CurveWriter>,
    pub writer_r: ManuallyDrop<CurveWriter>,
    pub writer_i: ManuallyDrop<CurveWriter>,
    pub writer_ever: ManuallyDrop<CurveWriter>,
    pub paths: [String; 4]
}

impl Drop for SirWriter
{
    fn drop(&mut self)
    {
        // first drop all the Writers!
        unsafe{
            ManuallyDrop::drop(&mut self.writer_s);
            ManuallyDrop::drop(&mut self.writer_r);
            ManuallyDrop::drop(&mut self.writer_i);
            ManuallyDrop::drop(&mut self.writer_ever);
        };

        // next: Zipping time!
        for path in self.paths.iter()
        {
            let out = Command::new("gzip")
                .arg(path)
                .output();
            match out {
                Ok(_) => println!("Success! Zipped {path}"),
                Err(e) => println!("Error! Failed to zip {path} due to {:?}", e)
            }
        }
    }
}

impl SirWriter
{
    #[inline]
    pub fn writer_iter(&mut self) -> impl Iterator<Item=&mut CurveWriter>
    {
        let slice = [
            self.writer_s.deref_mut(), 
            self.writer_r.deref_mut(), 
            self.writer_i.deref_mut(), 
            self.writer_ever.deref_mut()
        ];
        slice.into_iter()
    }

    pub fn new(name: &str, number: usize) -> Self
    {
        let names: [String; 4] = [
            format!("{name}_{number}_s.curves"),
            format!("{name}_{number}_r.curves"),
            format!("{name}_{number}_i.curves"),
            format!("{name}_{number}_e.curves")
        ];

        let mut files = names.clone().map(
            |name| 
            {
                BufWriter::new(
                    File::create(name)
                        .expect("unable to create file S")
                )
            }
        ).into_iter();


        Self{
            writer_s: ManuallyDrop::new(files.next().unwrap()),
            writer_r: ManuallyDrop::new(files.next().unwrap()),
            writer_i: ManuallyDrop::new(files.next().unwrap()),
            writer_ever: ManuallyDrop::new(files.next().unwrap()),
            paths: names
        }
    }

    pub fn write_energy<E, I>(&mut self, energy: E, extinction_index: I) -> std::io::Result<()>
    where E: Display,
        I: Display
    {
        for w in self.writer_iter()
        {
            write!(w, "{energy} {extinction_index} ")?;
        }
        
        Ok(())
    }

    pub fn write_current<A>(&mut self, graph: &GenericGraph<InfectionState, A>) -> std::io::Result<()>
    where A: AdjContainer<InfectionState>
    {
        let mut i = 0_u32;
        let mut r = 0_u32;
        let mut s = 0_u32;
        
        graph.contained_iter()
            .for_each(
                |contained|
                match contained{
                    InfectionState::Infected => i +=1,
                    InfectionState::Recovered => r +=1,
                    InfectionState::Suspectible => s += 1
                    
                }
            );

        write!(self.writer_i, "{i} ")?;
        write!(self.writer_r, "{r} ")?;
        write!(self.writer_s, "{s} ")?;
        let e = i + r;
        write!(self.writer_ever, "{e} ")
    }

    pub fn write_line(&mut self) -> std::io::Result<()>
    {
        for w in self.writer_iter()
        {
            writeln!(w)?;
        }
        Ok(())
    }

    pub fn write_header(&mut self, jsons: &[Value]) -> std::io::Result<()>
    {
        for w in self.writer_iter()
        {
            write_header(w)?;
            write_jsons(jsons, w)?;
        }

        Ok(())
    }
}

pub fn write_jsons<W: Write>(jsons: &[Value], mut writer: W) -> std::io::Result<()>
{
    for j in jsons{
        write!(writer, "#")?;
        serde_json::to_writer(&mut writer, j)?;
        writeln!(writer)?;
    }
    Ok(())
}

fn write_header(writer: &mut CurveWriter) -> std::io::Result<()>
{
    writeln!(writer, "#Energy Extinction_index Curve[0] Curve[1] …")?;
    Ok(())
}