% F11*4,F9(starting 'OMR pre-AHC vs AHC')
names = {'loomstart','RB_reg0.5_inRB','OMR (OT=non_preTec + preTec)','PT','OMR pre-AHC vs AHC','RB_reg0.5'      ...
    };
OT_left = [4496,219,62+61,273,252,384         ...
    ];
OT_total = [7505,526,131+297,591,502,850      ...
    ];

AH_left = [1112,262,85,184,292,193 ...
    ];
AH_total = [1409,1992,528,651,763,1288       ...
    ];
    
% preTec_left = [61,
%     ];
% preTec_total = [287 ...
%     ];
%%
OT_left./OT_total
AH_left./AH_total

abs((AH_left./AH_total-0.5)./(OT_left./OT_total - 0.5))