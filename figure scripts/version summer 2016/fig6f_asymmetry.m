
% OTec_left./OTec_total:   0.5991    0.4163    0.2874    0.4619    0.5020    0.4518
% AH_left./AH_total:          0.7892    0.1315    0.1610    0.2826    0.3827    0.1498
% abs((AH_left./AH_total-0.5)./(OTec_left./OTec_total - 0.5)): 2.9194    4.4049    1.5945    5.7093   58.8847    7.2593


 m1 = [0.5991,0.4163,0.2874,0.4619,0.5020,0.4518];
 m2 = [0.7892,0.1315,0.1610,0.2826,0.3827,0.1498];
 m3 = [2.9194,4.4049,1.5945,5.7093,58.8847,7.2593];
 y = vertcat(m1,m2)';
 y2 = vertcat(abs(m1-0.5)+0.5,abs(m2-0.5)+0.5)';
 
 figure;
 bar(y,'BaseValue',0.5);xlim([0,length(m1)+1]);ylim([0,1]);
 ylabel('left/total')
 xlabel('fish #')
 legend('upstream','downstream')
%  bar(y2,'BaseValue',0.5);ylim([0.5,1])


%%
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