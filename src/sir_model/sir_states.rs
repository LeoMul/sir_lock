use {
    serde::{Serialize, Deserialize},
    net_ensembles::Node
};

#[derive(Clone, Debug,PartialEq,Copy)]
#[derive(Serialize, Deserialize)]
pub enum InfectionState{
    Suspectible,
    Infected,
    Recovered,
}
impl InfectionState{
    pub fn sus_check(&self) -> bool{
        matches!(self,InfectionState::Suspectible)
    }
    pub fn inf_check(&self) -> bool{
        matches!(self,InfectionState::Infected)
    }
    pub fn rec_check(&self) -> bool{
        matches!(self,InfectionState::Recovered)
    }
    
    pub fn is_or_was_infected(&self) -> bool
    {
        matches!(self, Self::Infected | Self::Recovered)
    }
    }
    



impl Default for InfectionState{
    fn default() -> Self{
        InfectionState::Suspectible
    }
}
impl Node for InfectionState{
    fn new_from_index(_index: usize) -> Self{
        InfectionState::Suspectible
    }

}