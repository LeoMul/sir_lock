self.base_model.ensemble.contained_iter()
            .zip(
                self.lock_graph.contained_iter_mut()
            ).for_each(
                |(old,new)|
                {
                    *new=*old
                }
            );