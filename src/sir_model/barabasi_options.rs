

use{
    std::num::*,
    crate::scan_lambda::*,
    crate::scan_gamma::*,
    crate::scan_lambda_gamma::*,
    crate::GraphType,
    //crate::scan_lambda_gamma::*,
    //crate::time_graph::*,
    //crate::lifespanhist::*,
    crate::life_span_size_fitting_lambdascan::*,
    crate::life_span_size_fitting_threshscan::*,

    crate::prop_test::*,
    crate::scan_lambda_lock_thresh::*,
    crate::critical_lambda::*,
    crate::simple_sampling::*,
    crate::scan_lock_params::*,
    crate::critical_threshold::*,
    crate::connectedcomponent::*,
    crate::simplecurves::*
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
    pub fn from_connectedcomponent_param(param: &ConnectedComponentParams,seed:u64) -> Self
    {
        let (m,source_n )= match param.graph_type {
            GraphType::Barabasi(mm,source_nn) => (mm,source_nn),
            _ => panic!("Invalid graph type")
        };
        Self{
            m,
            source_n,
            graph_seed: seed,
            system_size: param.system_size,
            lambda: 0.,
            gamma: 0.
        }
    }
    pub fn from_lock_scan_param(param: &ScanLockParams) -> Self
    {
        let (m,source_n )= match param.graph_type {
            GraphType::Barabasi(mm,source_nn) => (mm,source_nn),
            _ => panic!("Invalid graph type")
        };
        Self{
            m,
            source_n,
            graph_seed: param.graph_seed,
            system_size: param.system_size,
            lambda: param.trans_prob,
            gamma: param.recovery_prob
        }
    }
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
    pub fn from_gamma_scan_param(param: &ScanGammaParams) -> Self{
        let (m,source_n )= match param.graph_type {
            GraphType::Barabasi(mm,source_nn) => (mm,source_nn),
            _ => panic!("Invalid graph type")
        };
        Self{
            graph_seed: param.graph_seed,
            system_size: param.system_size,
            lambda: param.transmission_prob,
            gamma: param.gamma_range.start,
            m,
            source_n,
        }
    }
    pub fn from_lambda_gamma_scan_param(param: &ScanLambdaGammaParams) -> Self{
        let (m,source_n )= match param.graph_type {
            GraphType::Barabasi(mm,source_nn) => (mm,source_nn),
            _ => panic!("Invalid graph type")
        };
        Self{
            graph_seed: param.graph_seed,
            system_size: param.system_size,
            lambda: param.lambda_range.start,
            gamma: param.gamma_range.start,
            m,
            source_n,
        }

    }
    pub fn from_lambda_thresh_param(param: &ScanLambdaThreshParams) -> Self{
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
    pub fn from_critical_lambda_params(param: &CriticalLambdaParams,system_size_new:NonZeroUsize,graph_seed:u64) -> Self{
        let (m,source_n )= match param.graph_type {
            GraphType::Barabasi(mm,source_nn) => (mm,source_nn),
            _ => panic!("Invalid graph type")
        };
        Self{
            graph_seed,
            system_size: system_size_new,
            lambda: param.lambda_range.start,
            gamma: param.recovery_prob,
            m,
            source_n,
        }

    }
    pub fn from_critical_thresh_params(param: &CriticalThreshParams,system_size_new:NonZeroUsize,graph_seed:u64) -> Self{
        let (m,source_n )= match param.graph_type {
            GraphType::Barabasi(mm,source_nn) => (mm,source_nn),
            _ => panic!("Invalid graph type")
        };
        Self{
            graph_seed,
            system_size: system_size_new,
            lambda: param.lambda,
            gamma: param.recovery_prob,
            m,
            source_n,
        }

    }
    pub fn from_life_span_fiting_params(param: &LifespanSizeFittingParams,system_size_new:NonZeroUsize,graph_seed:u64) -> Self{
        let (m,source_n )= match param.graph_type {
            GraphType::Barabasi(mm,source_nn) => (mm,source_nn),
            _ => panic!("Invalid graph type")
        };
        Self{
            graph_seed,
            system_size: system_size_new,
            lambda: param.trans_prob_range.start,
            gamma: param.recovery_prob,
            m,
            source_n,
        }

    }
    pub fn from_life_span_fiting_thresh_params(param: &LifespanSizeFittingThreshParams,system_size_new:NonZeroUsize,graph_seed:u64) -> Self{
        let (m,source_n )= match param.graph_type {
            GraphType::Barabasi(mm,source_nn) => (mm,source_nn),
            _ => panic!("Invalid graph type")
        };
        Self{
            graph_seed,
            system_size: system_size_new,
            lambda: param.trans_prob,
            gamma: param.recovery_prob,
            m,
            source_n,
        }

    }
    pub fn from_simple_sample(param: &SimpleSampleParam) -> Self
    {
        let (m,source_n) = match param.graph_type {
            GraphType::Barabasi(mm,source_nn) => (mm,source_nn),
            _ => panic!("Invalid graph type")
        };
        Self{
            graph_seed: param.graph_seed,
            system_size: param.system_size,
            lambda: param.lambda,
            gamma: param.recovery_prob,
            m,
            source_n
        }
    }
    pub fn from_simplecurves(param: &SimpleCurvesParam) -> Self
    {
        let (m,source_n) = match param.graph_type {
            GraphType::Barabasi(mm,source_nn) => (mm,source_nn),
            _ => panic!("Invalid graph type")
        };
        Self{
            graph_seed: param.graph_seed,
            system_size: param.system_size,
            lambda: param.lambda,
            gamma: param.recovery_prob,
            m,
            source_n
        }
    }


}