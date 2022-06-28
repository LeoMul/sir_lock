use std::io::Write;
use super::*;


#[derive(Clone)]
pub struct GridMapF64Generic<T>
{
    grid: GridF64,
    values: Vec<T>
}

impl <T> GridMapF64Generic<T>{
    pub fn from_fn<F>(grid: GridF64, mapper: F) -> Self
    where F: FnMut (Point2D) -> T
    {
        let vec: Vec<_> = grid.grid_point2d_iter()
            .map(mapper)
            .collect();
        Self{
            grid,
            values: vec
        }
    }
    pub fn from_vec_unchecked(grid: GridF64, z: Vec<T>) -> Self
    {
        Self{
            grid,
            values: z
        }
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

}