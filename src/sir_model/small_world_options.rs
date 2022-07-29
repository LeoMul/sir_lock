use{
    std::num::*,
    crate::scan_lambda::*,
    crate::scan_gamma::*,
    crate::GraphType,
    crate::scan_lambda_gamma::*,
    //crate::time_graph::*,
    crate::lifespanhist::*,
    crate::life_span_size_fitting::*,
    crate::critical_lambda::*,
    crate::simple_sampling::*,
    crate::prop_test::*,
    crate::scan_lambda_lock_thresh::*,
    crate::scan_lock_params::*,
    crate::critical_threshold::*
};

pub struct SWOptions{
    pub graph_seed: u64,
    pub rewire_prob: f64,
    pub system_size: NonZeroUsize,
    pub lambda: f64,
    pub gamma: f64
}

impl SWOptions{
    pub fn from_lock_scan_param(param: &ScanLockParams) -> Self
    {
        let rewire_prob = match param.graph_type {
            GraphType::SmallWorld(rewire) => rewire,
            _ => panic!("Invalid graph type")
        };
        Self{
            rewire_prob,
            graph_seed: param.graph_seed,
            system_size: param.system_size,
            lambda: param.trans_prob,
            gamma: param.recovery_prob
        }
    }
    pub fn from_lambda_scan_param(param: &ScanLambdaParams) -> Self
    {
        let rewire_prob = match param.graph_type {
            GraphType::SmallWorld(rewire) => rewire,
            _ => panic!("Invalid graph type")
        };
        Self{
            rewire_prob,
            graph_seed: param.graph_seed,
            system_size: param.system_size,
            lambda: param.lambda_range.start,
            gamma: param.recovery_prob
        }
    }
    pub fn from_gamma_scan_param(param: &ScanGammaParams) -> Self
    {
        let rewire_prob = match param.graph_type {
            GraphType::SmallWorld(rewire) => rewire,
            _ => panic!("Invalid graph type")
        };
        Self{
            rewire_prob,
            graph_seed: param.graph_seed,
            system_size: param.system_size,
            lambda: param.transmission_prob,
            gamma: param.gamma_range.start
        }
    }
    pub fn from_lambda_gamma_scan_param(param: &ScanLambdaGammaParams) -> Self
    {
        let rewire_prob = match param.graph_type {
            GraphType::SmallWorld(rewire) => rewire,
            _ => panic!("Invalid graph type")
        };
        Self{
            rewire_prob,
            graph_seed: param.graph_seed,
            system_size: param.system_size,
            lambda: param.lambda_range.start,
            gamma: param.gamma_range.start,
        }
    }
   // pub fn from_time_graph_param(param: &TimeGraphParams) -> Self{
     //   let rewire_prob = match param.graph_type {
       //     GraphType::SmallWorld(rewire) => rewire,
         //   _ => panic!("Invalid graph type")
        //};
        //Self{
          //  rewire_prob,
            //graph_seed: param.graph_seed,
            //system_size: param.system_size,
            //lambda: param.trans_prob,
        //    gamma: param.recovery_prob
        //}

    //}
    pub fn from_lifespan_param(param: &LifeSpanParams) -> Self{
        let rewire_prob = match param.graph_type {
            GraphType::SmallWorld(rewire) => rewire,
            _ => panic!("Invalid graph type")
        };
        Self{
            rewire_prob,
            graph_seed: param.graph_seed,
            system_size: param.system_size,
            lambda: param.trans_prob,
            gamma: param.recovery_prob
        }

    }
    pub fn from_lifespan_size_fitting_param(param: &LifespanSizeFittingParams,system_size_new:NonZeroUsize,new_graph_seed:u64) -> Self{
        let rewire_prob = match param.graph_type {
            GraphType::SmallWorld(rewire) => rewire,
            _ => panic!("Invalid graph type")
        };
        Self{
            rewire_prob,
            graph_seed: new_graph_seed,
            system_size: system_size_new,
            lambda: param.trans_prob_range.start,
            gamma: param.recovery_prob
        }

    }
    pub fn from_lifespan_size_fitting_param_and_size(param: &LifespanSizeFittingParams,system_size_new:NonZeroUsize) -> Self{
        let rewire_prob = match param.graph_type {
            GraphType::SmallWorld(rewire) => rewire,
            _ => panic!("Invalid graph type")
        };
        Self{
            rewire_prob,
            graph_seed: param.graph_seed,
            system_size: system_size_new,
            lambda: param.trans_prob_range.start,
            gamma: param.recovery_prob
        }

    }

    pub fn from_critical_lambda_params(param:&CriticalLambdaParams,system_size_new:NonZeroUsize,new_graph_seed:u64) -> Self{
        let rewire_prob = match param.graph_type {
            GraphType::SmallWorld(rewire) => rewire,
            _ => panic!("Invalid graph type")
        };
        Self{
            rewire_prob,
            graph_seed: new_graph_seed,
            system_size: system_size_new,
            lambda: param.lambda_range.start,
            gamma: param.recovery_prob
        }

    }
    pub fn from_critical_thresh_params(param:&CriticalThreshParams,system_size_new:NonZeroUsize,new_graph_seed:u64) -> Self{
        let rewire_prob = match param.graph_type {
            GraphType::SmallWorld(rewire) => rewire,
            _ => panic!("Invalid graph type")
        };
        Self{
            rewire_prob,
            graph_seed: new_graph_seed,
            system_size: system_size_new,
            lambda: param.lambda,
            gamma: param.recovery_prob
        }

    }
    pub fn from_simple_sample(param: &SimpleSampleParam) -> Self
    {
        let rewire_prob = match param.graph_type {
            GraphType::SmallWorld(rewire) => rewire,
            _ => panic!("Invalid graph type")
        };
        Self{
            graph_seed: param.graph_seed,
            system_size: param.system_size,
            lambda: param.lambda,
            gamma: param.recovery_prob,
            rewire_prob,
        }
    }
    pub fn from_prop_test_param(param:&PropTestParams) -> Self{
        let rewire_prob = match param.graph_type {
            GraphType::SmallWorld(rewire) => rewire,
            _ => panic!("Invalid graph type")
        };
        Self{
            graph_seed: param.graph_seed,
            system_size: param.system_size,
            lambda: param.lambda,
            gamma: param.recovery_prob,
            rewire_prob,}

    }
    pub fn from_lambda_thresh_param(param: &ScanLambdaThreshParams) -> Self{
        let rewire_prob = match param.graph_type {
            GraphType::SmallWorld(rewire) => rewire,
            _ => panic!("Invalid graph type")
        };
        Self{
            graph_seed: param.graph_seed,
            system_size: param.system_size,
            lambda: param.lambda_range.start,
            gamma: param.recovery_prob,
            rewire_prob
        }
    }
    
}