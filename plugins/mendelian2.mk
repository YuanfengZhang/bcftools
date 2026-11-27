# Track the policy headers so incremental builds cannot retain an old evaluator.
plugins/mendelian2$(PLUGIN_EXT): plugins/mendelian2-xy.h plugins/mendelian2-xy-eval.h
