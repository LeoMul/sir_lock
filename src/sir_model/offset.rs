use{
    serde::{Serialize, Deserialize},
    std::num::*
};



/// # Index Mapping for RNG Vectors
/// * This does the `rotation` of the rng vecs without actually 
/// copying the whole (and very large) list - this is much more efficient
#[derive(Serialize, Deserialize, Debug, Clone, Copy)]
pub struct Offset{
    offset: usize,
    n: usize,
    bound: NonZeroUsize,
    time_with_offset: usize,
}

impl Offset {

    /// # Set the current time
    /// * has to be called before `self.lookup_index`
    /// * This is an optimization - It was once part of the `lookup_index`
    ///  function - the way it is now has to be calculated only once per timestep, not for 
    /// each node
    #[inline(always)]
    pub fn set_time(&mut self, time: usize) {
        self.time_with_offset = time + self.offset;
        if self.time_with_offset >= self.bound.get() {
            self.time_with_offset -= self.bound.get();
        }

        // I believe, if a panic happens, which would be fixed by the following, then
        // something else is wrong in the programm and therefore, the following fix
        // should never be used!
        // 
        // self.time_with_offset = self.time_with_offset % self.bound;
    }

    /// Increase offset by 1 - i.e., Rotate RNG Vectors 
    /// * Wraps around bounds
    pub fn plus_1(&mut self){
        self.offset += 1;
        if self.offset >= self.bound.get() {
            self.offset = 0;
        }
    }

    /// Decrease offset by 1 - i.e., Rotate RNG Vectors in oposite direction
    pub fn minus_1(&mut self){
        self.offset = self.offset.checked_sub(1)
            .unwrap_or(self.bound.get() - 1);
    }

    pub fn new(bound: usize, n: usize) -> Self{
        Self{
            offset: 0,
            n,
            bound: unsafe{NonZeroUsize::new_unchecked(bound.max(1))},
            time_with_offset: 0
        }
    }

    /// choose a new bound
    pub fn set_bound(&mut self, bound: usize)
    {
        self.bound = unsafe{NonZeroUsize::new_unchecked(bound.max(1))};
    }

    /// You have to call `set_time` to set the correct time.
    /// If the offset changed, `set_time` has to be called again,
    /// to guarantee correct results.
    #[inline]
    pub fn lookup_index(&self, index: usize) -> usize
    {
        self.time_with_offset * self.n + index
    }

    pub fn get_time_with_offset(&self) -> usize
    {
        self.time_with_offset
    }
}