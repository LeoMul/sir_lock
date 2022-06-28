

use{
    std::num::*,
    crate::scan_lambda::*,
    crate::GraphType,
    //crate::scan_lambda_gamma::*,
    //crate::time_graph::*,
    //crate::lifespanhist::*,
    //crate::life_span_size_fitting::*,
   // crate::critical_lambda::*,
    crate::prop_test::*
};

pub struct BarabasiOptions{
    pub graph_seed: u64,
    pub system_size: NonZeroUsize,
    pub m:usize,
    pub source_n:usize,
    pub lambda: f64,
    pub gamma: f64,
   
    
}
impl BarabasiOptions{
    pub fn from_lambda_scan_param(param: &ScanLambdaParams) -> Self{
        let (m,source_n )= match param.graph_type {
            GraphType::Barabasi(mm,source_nn) => (mm,source_nn),
            _ => panic!("Invalid graph type")
        };
        Self{
            graph_seed: param.graph_seed,
            system_size: param.system_size,
            lambda: param.lambda_range.start,
            gamma: param.recovery_prob,
            m,
            source_n,
        }
    }
    pub fn from_prop_test_param(param:&PropTestParams) -> Self{
        let (m,source_n )= match param.graph_type {
            GraphType::Barabasi(mm,source_nn) => (mm,source_nn),
            _ => panic!("Invalid graph type")
        };
        Self{
            graph_seed: param.graph_seed,
            system_size: param.system_size,
            lambda: param.lambda,
            gamma: param.recovery_prob,
            m,
            source_n,
        }

    }
}