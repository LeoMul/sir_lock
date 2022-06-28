use std::{ io::Write};
use{
    serde::{Serialize,  de::DeserializeOwned},
    serde_json::Value,
    std::{
        
        process::exit,
        fs::File,
        io::BufReader,
        
        path::Path
    },
    
    
};


//Credit: Yannick Feld, this function originally written by him.

pub fn parse<P, T>(file: Option<P>) -> (T, Value)
where P: AsRef<Path>,
    T: Default + Serialize + DeserializeOwned
{
    match file
    {
        None => {
            let example = T::default();
            serde_json::to_writer_pretty(
                std::io::stdout(),
                &example
            ).expect("Unable to reach stdout");
            exit(0)
        }, 
        Some(file) => {
            let f = File::open(file)
                .expect("Unable to open file");
            let buf = BufReader::new(f);

            let json_val: Value = match serde_json::from_reader(buf)
            {
                Ok(v) => v,
                Err(e) => {
                    eprintln!("json parsing error!");
                    dbg!(e);
                    exit(1);
                }
            };

            let opt: T = match serde_json::from_value(json_val.clone())
            {
                Ok(o) => o,
                Err(e) => {
                    eprintln!("json parsing error!");
                    dbg!(e);
                    exit(1);
                }
            };

            (opt, json_val)    
        }
    }
}

pub fn write_json<W: Write>(mut writer: W, json: &Value)
{
    write!(writer, "#").unwrap();
    serde_json::to_writer(&mut writer, json).unwrap();
    writeln!(writer).unwrap();
}