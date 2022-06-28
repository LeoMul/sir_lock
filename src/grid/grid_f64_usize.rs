use std::{ops::RangeInclusive, num::NonZeroUsize};

use super::*;

#[derive(Debug, Clone, Copy)]
pub struct Point2D{
    pub x: f64,
    pub y: f64
}

#[derive(Clone)]
pub struct GridF64Usize{
    x_range: GridRangeF64,
    y_range: RangeInclusive<usize>,
    step_size_y: NonZeroUsize
}

impl GridF64Usize{
    pub fn new(
        x_range: GridRangeF64,
        y_range: RangeInclusive<usize>,
        step_size_y: NonZeroUsize
    ) -> Self
    {
        Self{
            x_range,
            y_range,
            step_size_y
        }
    }

    pub fn x_range(&self) -> &GridRangeF64
    {
        &self.x_range
    }

    pub fn x_range_iter(&self) -> GridRangeIterF64
    {
        self.x_range.iter()
    }

    pub fn y_range(&self) -> &RangeInclusive<usize>
    {
        &self.y_range
    }

    pub fn y_range_iter(&self) -> impl Iterator<Item=usize>
    {
        self.y_range
            .clone()
            .step_by(self.step_size_y.get())
    }


    pub fn contains_xy(&self, x: f64, y: usize) -> bool
    {
        self.x_range.contains(&x) && self.y_range.contains(&y)
    }

    pub fn grid_iter(&'_ self) -> impl Iterator<Item=(f64, usize)> + '_
    {
        self.x_range.iter()
            .flat_map(
                |x| 
                {
                    std::iter::repeat(x)
                        .zip(self.y_range_iter())
                }
            )
    }

    pub fn grid_point2d_iter(&'_ self) -> impl Iterator<Item=Point2DGeneric<f64, usize>> + '_
    {
        self.grid_iter()
            .map(|(x, y)| Point2DGeneric::<f64, usize>{x, y})
    }
}

#[cfg(test)]
mod testing
{
    use super::*;

    #[test]
    fn iter_test()
    {
        let x_range = GridRangeF64::new(0.0, 1.0, 3);
        let y_range = GridRangeF64::new(4.0, 3.0, 3);

        let grid = GridF64::new(x_range, y_range);

        let mut iter = grid.grid_iter();

        assert_eq!(Some((0.0, 4.0)), iter.next());
        assert_eq!(Some((0.0, 3.5)), iter.next());
        assert_eq!(Some((0.0, 3.0)), iter.next());
        assert_eq!(Some((0.5, 4.0)), iter.next());
        assert_eq!(Some((0.5, 3.5)), iter.next());
        assert_eq!(Some((0.5, 3.0)), iter.next());
        assert_eq!(Some((1.0, 4.0)), iter.next());
        assert_eq!(Some((1.0, 3.5)), iter.next());
        assert_eq!(Some((1.0, 3.0)), iter.next());
        assert_eq!(None, iter.next());
    }
}