use std::io::Write;
use super::*;

//pub type UsizeYPoint3D<T> = Point3DGeneric<f64, usize, T>;

#[derive(Clone)]
pub struct GridMapUsizeF64<T>
{
    grid: GridF64Usize,
    values: Vec<T>
}

impl<T> GridMapUsizeF64<T>{

    pub fn from_fn<F>(grid: GridF64Usize, mapper: F) -> Self
    where F: FnMut (Point2DGeneric<f64, usize>) -> T
    {
        let vec: Vec<_> = grid.grid_point2d_iter()
            .map(mapper)
            .collect();
        Self{
            grid,
            values: vec
        }
    }

    pub fn from_vec_unchecked(grid: GridF64Usize, z: Vec<T>) -> Self
    {
        Self{
            grid,
            values: z
        }
    }

    pub fn iter(&'_ self) -> impl Iterator<Item=((f64, usize), &T)> + '_
    {
        self.grid.grid_iter()
            .zip(self.values.iter())
    }


    pub fn write<W, W2>(&self, writer: &mut W, mut to_f64: W2) -> std::io::Result<()>
    where W: Write,
        W2: FnMut (&T) -> f64,
    {
        writeln!(writer, "#X Y Z")?;

        let iter = self.grid.x_range()
            .iter()
            .flat_map(
                |x| 
                {
                    std::iter::once(true)
                        .chain(std::iter::repeat(false))
                        .zip(std::iter::repeat(x))
                        .zip(self.grid.y_range_iter())
                }
            ).zip(self.values.iter());
        for (((new_line, x), y), z) in iter
        {
            if new_line {
                writeln!(writer)?;
            }
            writeln!(writer, "{:E} {:E} {:E}", x, y, to_f64(z))?;
        }
        Ok(())
    }

    /// # Create Something that can be plotted with Gnuplot!
    /// In the gnuplotterminal (opend by writing gnuplot in the Terminal) write `load filename` to load
    pub fn write_gnuplot<W, W2>(&
        self, 
        mut writer: W, 
        to_f64: W2,
        data_name: &str
    ) -> std::io::Result<()>
    where W: Write,
    W2: FnMut (&T) -> f64
    {
        writeln!(writer, "${data_name} << EOD")?;
        self.write(&mut writer, to_f64)?;
        writeln!(writer, "EOD")?;
        writeln!(writer, "splot ${data_name} with lines")
    }
}

#[cfg(test)]
mod grid_val_tests
{
    use super::*;
    use std::fs::File;
    use std::io::BufWriter;

    #[test]
    fn writer()
    {
        let file = File::create("grid_test.gp")
            .unwrap();
        let buf = BufWriter::new(file);

        let x_range = GridRangeF64::new(1.0,-1.0, 9);
        let y_range = GridRangeF64::new(-2.0, 2.0, 15);

        let grid = GridF64::new(x_range, y_range);
        let grid_map = GridMapF64::from_fn(grid, |point| point.x * point.x + point.y * point.y);
        grid_map.write_gnuplot(buf).unwrap();
    }
}