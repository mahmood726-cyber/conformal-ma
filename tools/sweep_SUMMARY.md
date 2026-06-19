# sweep_SUMMARY

## Observed pattern

For tau2=0.05, conformal is not higher than standard for any random-effects law in the y_new_by_re_dist table (normal standard 0.9508 vs conformal 0.9133; t3 standard 0.9543 vs conformal 0.9109; skew standard 0.9548 vs conformal 0.9158; mixture standard 0.9493 vs conformal 0.9105); the mean conformal/standard width ratios are normal 1.012, t3 1.099, skew 1.018, mixture 0.98. For tau2=0.30, conformal is not higher than standard for any random-effects law in the y_new_by_re_dist table (normal standard 0.9255 vs conformal 0.9195; t3 standard 0.9271 vs conformal 0.9199; skew standard 0.938 vs conformal 0.9307; mixture standard 0.9237 vs conformal 0.9158); the mean conformal/standard width ratios are normal 1.115, t3 1.376, skew 1.135, mixture 1.041.

## MAIN RUN (tau2=0.05)

Source: `sim_output/honest_coverage_main.json`

### headline

```json
"headline": {
    "estimand": "y_new (next study at median precision)",
    "correct_model_normal_none": {
      "n": 16000,
      "y_new": {
        "standard": 0.9674,
        "hksj": 0.964,
        "conformal": 0.9283
      }
    },
    "misspecified_any": {
      "n": 304000,
      "y_new": {
        "standard": 0.9515,
        "hksj": 0.9466,
        "conformal": 0.9118
      }
    }
  }
```

### y_new_by_re_dist

```json
"y_new_by_re_dist": {
    "normal": {
      "n": 80000,
      "standard": 0.9508,
      "hksj": 0.9459,
      "conformal": 0.9133,
      "mean_width_ratio_conf_std": 1.012
    },
    "t3": {
      "n": 80000,
      "standard": 0.9543,
      "hksj": 0.9498,
      "conformal": 0.9109,
      "mean_width_ratio_conf_std": 1.099
    },
    "skew": {
      "n": 80000,
      "standard": 0.9548,
      "hksj": 0.9502,
      "conformal": 0.9158,
      "mean_width_ratio_conf_std": 1.018
    },
    "mixture": {
      "n": 80000,
      "standard": 0.9493,
      "hksj": 0.944,
      "conformal": 0.9105,
      "mean_width_ratio_conf_std": 0.98
    }
  }
```

### theta_new_by_re_dist_k

```json
"theta_new_by_re_dist_k": {
    "normal|5": {
      "n": 20000,
      "standard": 0.9937,
      "hksj": 0.9876,
      "conformal": 0.8934,
      "mean_width_ratio_conf_std": 0.63
    },
    "normal|10": {
      "n": 20000,
      "standard": 0.982,
      "hksj": 0.9807,
      "conformal": 0.9624,
      "mean_width_ratio_conf_std": 1.013
    },
    "normal|15": {
      "n": 20000,
      "standard": 0.9797,
      "hksj": 0.979,
      "conformal": 0.9797,
      "mean_width_ratio_conf_std": 1.147
    },
    "normal|25": {
      "n": 20000,
      "standard": 0.9804,
      "hksj": 0.9804,
      "conformal": 0.9901,
      "mean_width_ratio_conf_std": 1.258
    },
    "t3|5": {
      "n": 20000,
      "standard": 0.9917,
      "hksj": 0.9866,
      "conformal": 0.9012,
      "mean_width_ratio_conf_std": 0.615
    },
    "t3|10": {
      "n": 20000,
      "standard": 0.9806,
      "hksj": 0.9796,
      "conformal": 0.9642,
      "mean_width_ratio_conf_std": 1.059
    },
    "t3|15": {
      "n": 20000,
      "standard": 0.978,
      "hksj": 0.9778,
      "conformal": 0.978,
      "mean_width_ratio_conf_std": 1.255
    },
    "t3|25": {
      "n": 20000,
      "standard": 0.9784,
      "hksj": 0.9784,
      "conformal": 0.9875,
      "mean_width_ratio_conf_std": 1.468
    },
    "skew|5": {
      "n": 20000,
      "standard": 0.9941,
      "hksj": 0.9879,
      "conformal": 0.8928,
      "mean_width_ratio_conf_std": 0.632
    },
    "skew|10": {
      "n": 20000,
      "standard": 0.9865,
      "hksj": 0.9848,
      "conformal": 0.9682,
      "mean_width_ratio_conf_std": 1.022
    },
    "skew|15": {
      "n": 20000,
      "standard": 0.9853,
      "hksj": 0.9848,
      "conformal": 0.9853,
      "mean_width_ratio_conf_std": 1.151
    },
    "skew|25": {
      "n": 20000,
      "standard": 0.9892,
      "hksj": 0.9892,
      "conformal": 0.995,
      "mean_width_ratio_conf_std": 1.265
    },
    "mixture|5": {
      "n": 20000,
      "standard": 0.9929,
      "hksj": 0.9861,
      "conformal": 0.8835,
      "mean_width_ratio_conf_std": 0.624
    },
    "mixture|10": {
      "n": 20000,
      "standard": 0.9826,
      "hksj": 0.9812,
      "conformal": 0.961,
      "mean_width_ratio_conf_std": 0.995
    },
    "mixture|15": {
      "n": 20000,
      "standard": 0.9807,
      "hksj": 0.98,
      "conformal": 0.9797,
      "mean_width_ratio_conf_std": 1.107
    },
    "mixture|25": {
      "n": 20000,
      "standard": 0.9832,
      "hksj": 0.9831,
      "conformal": 0.9918,
      "mean_width_ratio_conf_std": 1.192
    }
  }
```

## SENSITIVITY (tau2=0.30)

Source: `sim_output/honest_coverage_tau2hi.json`

### headline

```json
"headline": {
    "estimand": "y_new (next study at median precision)",
    "correct_model_normal_none": {
      "n": 16000,
      "y_new": {
        "standard": 0.955,
        "hksj": 0.9494,
        "conformal": 0.9411
      }
    },
    "misspecified_any": {
      "n": 304000,
      "y_new": {
        "standard": 0.9272,
        "hksj": 0.9217,
        "conformal": 0.9204
      }
    }
  }
```

### y_new_by_re_dist

```json
"y_new_by_re_dist": {
    "normal": {
      "n": 80000,
      "standard": 0.9255,
      "hksj": 0.919,
      "conformal": 0.9195,
      "mean_width_ratio_conf_std": 1.115
    },
    "t3": {
      "n": 80000,
      "standard": 0.9271,
      "hksj": 0.922,
      "conformal": 0.9199,
      "mean_width_ratio_conf_std": 1.376
    },
    "skew": {
      "n": 80000,
      "standard": 0.938,
      "hksj": 0.9327,
      "conformal": 0.9307,
      "mean_width_ratio_conf_std": 1.135
    },
    "mixture": {
      "n": 80000,
      "standard": 0.9237,
      "hksj": 0.9188,
      "conformal": 0.9158,
      "mean_width_ratio_conf_std": 1.041
    }
  }
```

### theta_new_by_re_dist_k

```json
"theta_new_by_re_dist_k": {
    "normal|5": {
      "n": 20000,
      "standard": 0.962,
      "hksj": 0.9448,
      "conformal": 0.8736,
      "mean_width_ratio_conf_std": 0.832
    },
    "normal|10": {
      "n": 20000,
      "standard": 0.9334,
      "hksj": 0.9315,
      "conformal": 0.9312,
      "mean_width_ratio_conf_std": 1.127
    },
    "normal|15": {
      "n": 20000,
      "standard": 0.931,
      "hksj": 0.9308,
      "conformal": 0.9524,
      "mean_width_ratio_conf_std": 1.207
    },
    "normal|25": {
      "n": 20000,
      "standard": 0.9318,
      "hksj": 0.9318,
      "conformal": 0.9714,
      "mean_width_ratio_conf_std": 1.293
    },
    "t3|5": {
      "n": 20000,
      "standard": 0.9613,
      "hksj": 0.9456,
      "conformal": 0.8766,
      "mean_width_ratio_conf_std": 0.863
    },
    "t3|10": {
      "n": 20000,
      "standard": 0.9377,
      "hksj": 0.9358,
      "conformal": 0.9377,
      "mean_width_ratio_conf_std": 1.346
    },
    "t3|15": {
      "n": 20000,
      "standard": 0.9363,
      "hksj": 0.9362,
      "conformal": 0.958,
      "mean_width_ratio_conf_std": 1.537
    },
    "t3|25": {
      "n": 20000,
      "standard": 0.9359,
      "hksj": 0.9361,
      "conformal": 0.9709,
      "mean_width_ratio_conf_std": 1.757
    },
    "skew|5": {
      "n": 20000,
      "standard": 0.9683,
      "hksj": 0.9525,
      "conformal": 0.892,
      "mean_width_ratio_conf_std": 0.853
    },
    "skew|10": {
      "n": 20000,
      "standard": 0.9483,
      "hksj": 0.9463,
      "conformal": 0.9456,
      "mean_width_ratio_conf_std": 1.14
    },
    "skew|15": {
      "n": 20000,
      "standard": 0.9474,
      "hksj": 0.947,
      "conformal": 0.9628,
      "mean_width_ratio_conf_std": 1.228
    },
    "skew|25": {
      "n": 20000,
      "standard": 0.9477,
      "hksj": 0.9476,
      "conformal": 0.9783,
      "mean_width_ratio_conf_std": 1.321
    },
    "mixture|5": {
      "n": 20000,
      "standard": 0.946,
      "hksj": 0.9278,
      "conformal": 0.871,
      "mean_width_ratio_conf_std": 0.833
    },
    "mixture|10": {
      "n": 20000,
      "standard": 0.9339,
      "hksj": 0.9316,
      "conformal": 0.9338,
      "mean_width_ratio_conf_std": 1.058
    },
    "mixture|15": {
      "n": 20000,
      "standard": 0.9338,
      "hksj": 0.9333,
      "conformal": 0.955,
      "mean_width_ratio_conf_std": 1.106
    },
    "mixture|25": {
      "n": 20000,
      "standard": 0.9336,
      "hksj": 0.9336,
      "conformal": 0.9716,
      "mean_width_ratio_conf_std": 1.165
    }
  }
```
