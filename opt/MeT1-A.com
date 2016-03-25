g09 << +++ > MeT1-A.log
%mem=1900MB
%Nprocshared=4
%Chk=/short/d63/reimers/MeT1-A
#P ONIOM(cam-b3lyp:AMBER=hardfirst)=EmbedCharge iop(4/28=5)
   opt(z-matrix,maxcycle=7,nrscale) nosym scf(conver=6) 

  lysozyme partial opt MeT1-A

 -1 1  1 1  1 1
 N-N3-0.175 0 x0001 y0001 z0001    
 C-CT--0.024 0 x0002 y0002 z0002    
 C-C-0.597 0 x0003 y0003 z0003    
 O-O--0.568 0 x0004 y0004 z0004    
 C-CT-0.034 0 x0005 y0005 z0005    
 C-CT-0.002 0 x0006 y0006 z0006    
 S-S--0.274 0 x0007 y0007 z0007    
 C-CT--0.054 0 x0008 y0008 z0008    
 H-H-0.227 0 x0009 y0009 z0009    
 H-H1-0.088 0 x0010 y0010 z0010    
 H-HC-0.024 0 x0011 y0011 z0011    
 H-HC-0.024 0 x0012 y0012 z0012    
 H-H1-0.044 0 x0013 y0013 z0013    
 H-H1-0.044 0 x0014 y0014 z0014    
 H-H1-0.068 0 x0015 y0015 z0015    
 H-H1-0.068 0 x0016 y0016 z0016    
 H-H1-0.068 0 x0017 y0017 z0017    
 H-H-0.227 0 x0018 y0018 z0018    
 H-H-0.227 0 x0019 y0019 z0019    
 N-N--0.416 0 x0020 y0020 z0020    
 C-CT--0.087 0 x0021 y0021 z0021    
 C-C-0.597 0 x0022 y0022 z0022    
 O-O--0.568 0 x0023 y0023 z0023    
 C-CT-0.298 0 x0024 y0024 z0024    
 C-CT--0.319 0 x0025 y0025 z0025    
 C-CT--0.319 0 x0026 y0026 z0026    
 H-H-0.272 0 x0027 y0027 z0027    
 H-H1-0.097 0 x0028 y0028 z0028    
 H-HC--0.030 0 x0029 y0029 z0029    
 H-HC-0.079 0 x0030 y0030 z0030    
 H-HC-0.079 0 x0031 y0031 z0031    
 H-HC-0.079 0 x0032 y0032 z0032    
 H-HC-0.079 0 x0033 y0033 z0033    
 H-HC-0.079 0 x0034 y0034 z0034    
 H-HC-0.079 0 x0035 y0035 z0035    
 N-N--0.416 0 x0036 y0036 z0036    
 C-CT--0.025 0 x0037 y0037 z0037    
 C-C-0.597 0 x0038 y0038 z0038 M H-H1-0.029  
 O-O--0.568 0 x0039 y0039 z0039 M  
 H-H-0.272 0 x0040 y0040 z0040    
 H-H1-0.070 0 x0041 y0041 z0041    
 H-H1-0.070 0 x0042 y0042 z0042    
 N-N--0.416 0 x0043 y0043 z0043 M  
 C-CT--0.087 0 x0044 y0044 z0044 M  
 C-C-0.597 0 x0045 y0045 z0045 M  
 O-O--0.568 0 x0046 y0046 z0046 M  
 C-CT-0.298 0 x0047 y0047 z0047 M  
 C-CT--0.277 0 x0048 y0048 z0048 M  
 C-CT--0.319 0 x0049 y0049 z0049 M  
 H-H-0.272 0 x0050 y0050 z0050 M  
 H-H1-0.097 0 x0051 y0051 z0051 M  
 H-HC--0.030 0 x0052 y0052 z0052 M  
 H-HC-0.035 0 x0053 y0053 z0053 M  
 H-HC-0.073 0 x0054 y0054 z0054 M  
 H-HC-0.042 0 x0055 y0055 z0055 M  
 H-HC-0.079 0 x0056 y0056 z0056 M  
 H-HC-0.079 0 x0057 y0057 z0057 M  
 H-HC-0.079 0 x0058 y0058 z0058 M  
 N-N--0.416 0 x0059 y0059 z0059 M  
 C-CT--0.060 0 x0060 y0060 z0060 M  
 C-C-0.597 0 x0061 y0061 z0061 M  
 O-O--0.433 0 x0062 y0062 z0062 M  
 C-CT-0.130 0 x0063 y0063 z0063 M  
 C-CT--0.043 0 x0064 y0064 z0064 M  
 C-CT--0.320 0 x0065 y0065 z0065 M  
 C-CT--0.066 0 x0066 y0066 z0066 M  
 H-H-0.272 0 x0067 y0067 z0067 M  
 H-H1-0.087 0 x0068 y0068 z0068 M  
 H-HC-0.019 0 x0069 y0069 z0069 M  
 H-HC-0.024 0 x0070 y0070 z0070 M  
 H-HC-0.024 0 x0071 y0071 z0071 M  
 H-HC-0.088 0 x0072 y0072 z0072 M  
 H-HC-0.088 0 x0073 y0073 z0073 M  
 H-HC-0.088 0 x0074 y0074 z0074 M  
 H-HC-0.019 0 x0075 y0075 z0075 M  
 H-HC-0.019 0 x0076 y0076 z0076 M  
 H-HC-0.019 0 x0077 y0077 z0077 M  
 N-N--0.416 0 x0078 y0078 z0078 M  
 C-CT-0.034 0 x0079 y0079 z0079 M  
 C-C-0.412 0 x0080 y0080 z0080 M  
 O-O2--0.568 0 x0081 y0081 z0081 M  
 C-CT--0.182 0 x0082 y0082 z0082 M  
 O-O2--0.815 0 x0083 y0083 z0083 M  
 H-H-0.240 0 x0084 y0084 z0084 M  
 H-H1-0.082 0 x0085 y0085 z0085 M  
 H-HC-0.060 0 x0086 y0086 z0086 M  
 H-HC-0.060 0 x0087 y0087 z0087 M  
 H-HC-0.060 0 x0088 y0088 z0088 M  
 N-N3-0.175 0 x0089 y0089 z0089 M  
 C-CT--0.022 0 x0090 y0090 z0090 M  
 C-C-0.295 0 x0091 y0091 z0091 M  
 O-O--0.388 0 x0092 y0092 z0092 M  
 C-CT-0.034 0 x0093 y0093 z0093 M  
 C-CT-0.002 0 x0094 y0094 z0094 M  
 S-S--0.274 0 x0095 y0095 z0095 M  
 C-CT--0.054 0 x0096 y0096 z0096 M  
 H-H-0.227 0 x0097 y0097 z0097 M  
 H-H1-0.041 0 x0098 y0098 z0098 M  
 H-HC-0.024 0 x0099 y0099 z0099 M  
 H-HC-0.024 0 x0100 y0100 z0100 M  
 H-H1-0.044 0 x0101 y0101 z0101 M  
 H-H1-0.044 0 x0102 y0102 z0102 M  
 H-H1-0.068 0 x0103 y0103 z0103 M  
 H-H1-0.037 0 x0104 y0104 z0104 M  
 H-H1-0.068 0 x0105 y0105 z0105 M  
 H-H-0.227 0 x0106 y0106 z0106 M  
 H-H-0.185 0 x0107 y0107 z0107 M  
 N-N--0.416 0 x0108 y0108 z0108 M  
 C-CT--0.060 0 x0109 y0109 z0109 M  
 C-C-0.597 0 x0110 y0110 z0110 M  
 O-O--0.568 0 x0111 y0111 z0111 M  
 C-CT-0.130 0 x0112 y0112 z0112 M  
 C-CT--0.043 0 x0113 y0113 z0113 M  
 C-CT--0.320 0 x0114 y0114 z0114 M  
 C-CT--0.066 0 x0115 y0115 z0115 M  
 H-H-0.272 0 x0116 y0116 z0116 M  
 H-H1-0.087 0 x0117 y0117 z0117 M  
 H-HC-0.019 0 x0118 y0118 z0118 M  
 H-HC-0.024 0 x0119 y0119 z0119 M  
 H-HC-0.024 0 x0120 y0120 z0120 M  
 H-HC-0.088 0 x0121 y0121 z0121 M  
 H-HC-0.088 0 x0122 y0122 z0122 M  
 H-HC-0.088 0 x0123 y0123 z0123 M  
 H-HC-0.019 0 x0124 y0124 z0124 M  
 H-HC-0.019 0 x0125 y0125 z0125 M  
 H-HC-0.019 0 x0126 y0126 z0126 M  
 N-N--0.416 0 x0127 y0127 z0127 M  
 C-CT-0.034 0 x0128 y0128 z0128 M  
 C-C-0.412 0 x0129 y0129 z0129 M  
 O-O2--0.568 0 x0130 y0130 z0130 M  
 C-CT--0.182 0 x0131 y0131 z0131 M  
 O-O2--0.815 0 x0132 y0132 z0132 M  
 H-H-0.272 0 x0133 y0133 z0133 M  
 H-H1-0.082 0 x0134 y0134 z0134 M  
 H-HC-0.060 0 x0135 y0135 z0135 M  
 H-HC-0.060 0 x0136 y0136 z0136 M  
 H-HC-0.060 0 x0137 y0137 z0137 M  
 N-N--0.415 0 x0138 y0138 z0138 M  
 C-CT-0.034 0 x0139 y0139 z0139 M  
 C-C-0.412 0 x0140 y0140 z0140 M  
 O-O2--0.568 0 x0141 y0141 z0141 M  
 C-CT--0.182 0 x0142 y0142 z0142 M  
 O-O2--0.815 0 x0143 y0143 z0143 M  
 H-H-0.272 0 x0144 y0144 z0144 M  
 H-H1-0.082 0 x0145 y0145 z0145 M  
 H-HC-0.060 0 x0146 y0146 z0146 M  
 H-HC-0.060 0 x0147 y0147 z0147 M  
 H-HC-0.060 0 x0148 y0148 z0148 M  
 C-CT--0.025 0 x0149 y0149 z0149 M  
 H-H-0.272 0 x0150 y0150 z0150 M  
 N-N--0.416 0 x0151 y0151 z0151 M  
 H-C-0.597 0 x0152 y0152 z0152 M  
 H-H1-0.070 0 x0153 y0153 z0153 M  
 H-H1-0.070 0 x0154 y0154 z0154 M  
 C-CT--0.079 0 x0155 y0155 z0155 M  
 O-O--0.369 0 x0156 y0156 z0156 M  
 H-N--0.160 0 x0157 y0157 z0157 M  
 C-C-0.597 0 x0158 y0158 z0158 M  
 H-CT-0.205 0 x0159 y0159 z0159 M  
 H-H1-0.097 0 x0160 y0160 z0160 M  
 C-CT--0.010 0 x0161 y0161 z0161 M  
 H-H-0.038 0 x0162 y0162 z0162 M  
 N-N--0.077 0 x0163 y0163 z0163 M  
 H-C-0.063 0 x0164 y0164 z0164 M  
 H-CT-0.019 0 x0165 y0165 z0165 M  
 H-H1-0.014 0 x0166 y0166 z0166 M  
 C-CT--0.087 0 x0167 y0167 z0167 M  
 O-O--0.568 0 x0168 y0168 z0168 M  
 H-N--0.416 0 x0169 y0169 z0169 M  
 C-C-0.597 0 x0170 y0170 z0170 M  
 H-CT-0.298 0 x0171 y0171 z0171 M  
 H-H1-0.097 0 x0172 y0172 z0172 M  
 C-CT--0.009 0 x0173 y0173 z0173 M  
 O-O--0.165 0 x0174 y0174 z0174 M  
 H-N--0.055 0 x0175 y0175 z0175 M  
 C-C-0.234 0 x0176 y0176 z0176 M  
 H-CT-0.013 0 x0177 y0177 z0177 M  
 H-H1-0.015 0 x0178 y0178 z0178 M  
 Variables:
x0009 =     0.2653
y0009 =    14.7070
z0009 =     6.9492
x0010 =     1.4095
y0010 =    12.8757
z0010 =     8.2516
x0011 =    -0.8519
y0011 =    14.1299
z0011 =     9.8531
x0012 =     0.8190
y0012 =    14.7592
z0012 =     9.6533
x0013 =     1.7269
y0013 =    13.2168
z0013 =    11.1227
x0014 =     0.5117
y0014 =    11.9551
z0014 =    10.7240
x0015 =    -1.2875
y0015 =    15.6696
z0015 =    12.4698
x0016 =     0.0880
y0016 =    15.5890
z0016 =    11.3169
x0017 =     0.4010
y0017 =    15.6871
z0017 =    13.0833
x0018 =    -0.8704
y0018 =    13.7380
z0018 =     6.4940
x0019 =    -0.0681
y0019 =    13.2171
z0019 =     6.4125
 Constants:
x0001 =    -0.1630
y0001 =    13.8300
z0001 =     7.2090
x0002 =     0.3390
y0002 =    13.0740
z0002 =     8.4090
x0003 =    -0.4710
y0003 =    11.8510
z0003 =     8.5330
x0004 =    -1.6500
y0004 =    11.8510
z0004 =     8.2300
x0005 =     0.2030
y0005 =    13.8500
z0005 =     9.7160
x0006 =     0.6600
y0006 =    13.0240
z0006 =    10.9370
x0007 =    -0.2890
y0007 =    13.4710
z0007 =    12.4150
x0008 =    -0.2700
y0008 =    15.2830
z0008 =    12.3110
x0020 =     0.1850
y0020 =    10.8110
z0020 =     9.0070
x0021 =    -0.3470
y0021 =     9.4920
z0021 =     8.9950
x0022 =     0.3510
y0022 =     8.9510
z0022 =    10.2200
x0023 =     1.5840
y0023 =     8.9700
z0023 =    10.2830
x0024 =     0.1270
y0024 =     8.6880
z0024 =     7.7230
x0025 =    -0.4860
y0025 =     7.2310
z0025 =     7.6640
x0026 =    -0.1320
y0026 =     9.4110
z0026 =     6.4320
x0027 =     1.0985
y0027 =    11.0255
z0027 =     9.3807
x0028 =    -1.4454
y0028 =     9.4345
z0028 =     8.9821
x0029 =     1.2174
y0029 =     8.5995
z0029 =     7.8377
x0030 =    -0.6429
y0030 =     6.8581
z0030 =     8.6869
x0031 =     0.2076
y0031 =     6.5632
z0031 =     7.1321
x0032 =    -1.4484
y0032 =     7.2599
z0032 =     7.1321
x0033 =     0.8187
y0033 =     9.7860
z0033 =     6.0252
x0034 =    -0.8116
y0034 =    10.2562
z0034 =     6.6157
x0035 =    -0.5925
y0035 =     8.7195
z0035 =     5.7111
x0036 =    -0.4240
y0036 =     8.4720
z0036 =    11.1930
x0037 =     0.1590
y0037 =     7.9090
z0037 =    12.4230
x0038 =    -0.4980
y0038 =     6.6220
z0038 =    12.8910
x0039 =    -1.6670
y0039 =     6.3860
z0039 =    12.6440
x0040 =    -1.4146
y0040 =     8.5333
z0040 =    11.0055
x0041 =     0.0622
y0041 =     8.6566
z0041 =    13.2240
x0042 =     1.2239
y0042 =     7.7050
z0042 =    12.2378
x0043 =    -4.4140
y0043 =     7.1610
z0043 =    12.2700
x0044 =    -4.9890
y0044 =     8.3710
z0044 =    11.6340
x0045 =    -4.3360
y0045 =     8.9530
z0045 =    10.3610
x0046 =    -3.1100
y0046 =     8.8700
z0046 =    10.1430
x0047 =    -5.1490
y0047 =     9.6130
z0047 =    12.6310
x0048 =    -5.9980
y0048 =     9.2470
z0048 =    13.7550
x0049 =    -3.8070
y0049 =    10.1430
z0049 =    13.1310
x0050 =    -3.4338
y0050 =     6.9181
z0050 =    12.2851
x0051 =    -5.9438
y0051 =     7.9184
z0051 =    11.3282
x0052 =    -5.6217
y0052 =    10.4263
z0052 =    12.0609
x0053 =    -6.3459
y0053 =    10.1585
z0053 =    14.2631
x0054 =    -5.4242
y0054 =     8.6291
z0054 =    14.4613
x0055 =    -6.8648
y0055 =     8.6772
z0055 =    13.3889
x0056 =    -3.9153
y0056 =    11.1985
z0056 =    13.4211
x0057 =    -3.0581
y0057 =    10.0573
z0057 =    12.3299
x0058 =    -3.4823
y0058 =     9.5544
z0058 =    14.0017
x0059 =    -5.2020
y0059 =     9.5750
z0059 =     9.5490
x0060 =    -4.8330
y0060 =    10.5900
z0060 =     8.5540
x0061 =    -5.9950
y0061 =    11.5230
z0061 =     8.8230
x0062 =    -7.1170
y0062 =    11.0270
z0062 =     8.7910
x0063 =    -4.9830
y0063 =    10.0340
z0063 =     7.1170
x0064 =    -4.5840
y0064 =     8.5530
z0064 =     7.0580
x0065 =    -4.2140
y0065 =    10.9200
z0065 =     6.0910
x0066 =    -5.0540
y0066 =     7.7620
z0066 =     5.8560
x0067 =    -6.1550
y0067 =     9.2717
z0067 =     9.6902
x0068 =    -3.8098
y0068 =    10.9883
z0068 =     8.6210
x0069 =    -6.0434
y0069 =    10.0804
z0069 =     6.8283
x0070 =    -4.9896
y0070 =     8.0641
z0070 =     7.9560
x0071 =    -3.4853
y0071 =     8.5044
z0071 =     7.0781
x0072 =    -3.4565
y0072 =    11.5179
z0072 =     6.6189
x0073 =    -4.9225
y0073 =    11.5907
z0073 =     5.5829
x0074 =    -3.7210
y0074 =    10.2758
z0074 =     5.3481
x0075 =    -4.4736
y0075 =     6.8306
z0075 =     5.7809
x0076 =    -4.9082
y0076 =     8.3600
z0076 =     4.9443
x0077 =    -6.1212
y0077 =     7.5214
z0077 =     5.9705
x0078 =    -5.8360
y0078 =    12.8070
z0078 =     9.1770
x0079 =    -4.5880
y0079 =    13.4940
z0079 =     9.4650
x0080 =    -4.8320
y0080 =    14.3750
z0080 =    10.6680
x0081 =    -4.9360
y0081 =    15.5600
z0081 =    10.3570
x0082 =    -4.2280
y0082 =    14.4280
z0082 =     8.3030
x0083 =    -4.9450
y0083 =    14.0280
z0083 =    11.8620
x0084 =    -6.7204
y0084 =    13.2918
z0084 =     9.2306
x0085 =    -3.7858
y0085 =    12.7592
z0085 =     9.6280
x0086 =    -3.2380
y0086 =    14.1563
z0086 =     7.9080
x0087 =    -4.2075
y0087 =    15.4677
z0087 =     8.6617
x0088 =    -4.9807
y0088 =    14.3290
z0088 =     7.5070
x0089 =    -2.6230
y0089 =    18.3150
z0089 =     6.5530
x0090 =    -2.7640
y0090 =    18.9090
z0090 =     7.8970
x0091 =    -1.9800
y0091 =    20.1930
z0091 =     7.9460
x0092 =    -0.7510
y0092 =    20.2270
z0092 =     7.8210
x0093 =    -2.3320
y0093 =    17.9510
z0093 =     9.0300
x0094 =    -2.7660
y0094 =    18.4180
z0094 =    10.4610
x0095 =    -1.6240
y0095 =    19.4770
z0095 =    11.4670
x0096 =    -2.8430
y0096 =    20.2930
z0096 =    12.5320
x0097 =    -3.0847
y0097 =    17.4435
z0097 =     6.3351
x0098 =    -3.8315
y0098 =    19.1107
z0098 =     8.0695
x0099 =    -1.2354
y0099 =    17.8666
z0099 =     9.0118
x0100 =    -2.7794
y0100 =    16.9646
z0100 =     8.8379
x0101 =    -2.9580
y0101 =    17.5083
z0101 =    11.0488
x0102 =    -3.7025
y0102 =    18.9826
z0102 =    10.3419
x0103 =    -3.8046
y0103 =    19.7638
z0103 =    12.4600
x0104 =    -2.9727
y0104 =    21.3360
z0104 =    12.2076
x0105 =    -2.4903
y0105 =    20.2737
z0105 =    13.5737
x0106 =    -2.1573
y0106 =    18.5294
z0106 =     5.6828
x0107 =    -2.9671
y0107 =    18.9955
z0107 =     5.8915
x0108 =     4.2680
y0108 =     9.5750
z0108 =     9.5490
x0109 =     4.6370
y0109 =    10.5900
z0109 =     8.5540
x0110 =     3.4750
y0110 =    11.5230
z0110 =     8.8230
x0111 =     2.3530
y0111 =    11.0270
z0111 =     8.7910
x0112 =     4.4870
y0112 =    10.0340
z0112 =     7.1170
x0113 =     4.8860
y0113 =     8.5530
z0113 =     7.0580
x0114 =     5.2560
y0114 =    10.9200
z0114 =     6.0910
x0115 =     4.4160
y0115 =     7.7620
z0115 =     5.8560
x0116 =     3.3150
y0116 =     9.2717
z0116 =     9.6902
x0117 =     5.6602
y0117 =    10.9883
z0117 =     8.6210
x0118 =     3.4266
y0118 =    10.0804
z0118 =     6.8283
x0119 =     4.4804
y0119 =     8.0641
z0119 =     7.9560
x0120 =     5.9847
y0120 =     8.5044
z0120 =     7.0781
x0121 =     6.0135
y0121 =    11.5179
z0121 =     6.6189
x0122 =     4.5475
y0122 =    11.5907
z0122 =     5.5829
x0123 =     5.7490
y0123 =    10.2758
z0123 =     5.3481
x0124 =     4.9964
y0124 =     6.8306
z0124 =     5.7809
x0125 =     4.5618
y0125 =     8.3600
z0125 =     4.9443
x0126 =     3.3488
y0126 =     7.5214
z0126 =     5.9705
x0127 =     3.6340
y0127 =    12.8070
z0127 =     9.1770
x0128 =     4.8820
y0128 =    13.4940
z0128 =     9.4650
x0129 =     4.6380
y0129 =    14.3750
z0129 =    10.6680
x0130 =     4.5340
y0130 =    15.5600
z0130 =    10.3570
x0131 =     5.2420
y0131 =    14.4280
z0131 =     8.3030
x0132 =     4.5250
y0132 =    14.0280
z0132 =    11.8620
x0133 =     2.7496
y0133 =    13.2918
z0133 =     9.2306
x0134 =     5.6842
y0134 =    12.7592
z0134 =     9.6280
x0135 =     6.2320
y0135 =    14.1563
z0135 =     7.9080
x0136 =     5.2625
y0136 =    15.4677
z0136 =     8.6617
x0137 =     4.4893
y0137 =    14.3290
z0137 =     7.5070
x0138 =     2.0010
y0138 =    19.3720
z0138 =     7.7530
x0139 =     2.4880
y0139 =    18.0500
z0139 =     7.9910
x0140 =     1.9420
y0140 =    17.0480
z0140 =     7.0010
x0141 =     2.6950
y0141 =    16.1170
z0141 =     6.6810
x0142 =     2.1570
y0142 =    17.6330
z0142 =     9.4030
x0143 =     0.7860
y0143 =    17.1200
z0143 =     6.5320
x0144 =     1.0373
y0144 =    19.6183
z0144 =     7.5776
x0145 =     3.5796
y0145 =    18.0670
z0145 =     7.8562
x0146 =     1.0665
y0146 =    17.5372
z0146 =     9.5104
x0147 =     2.5299
y0147 =    18.3925
z0147 =    10.1060
x0148 =     2.6338
y0148 =    16.6658
z0148 =     9.6204
x0149 =    -0.3900
y0149 =     4.5590
z0149 =    14.0790
x0150 =     1.2012
y0150 =     6.0435
z0150 =    13.7266
x0151 =     0.2370
y0151 =     5.7780
z0151 =    13.5850
x0152 =     0.1271
y0152 =     4.1507
z0152 =    14.9598
x0153 =    -0.3749
y0153 =     3.8044
z0153 =    13.2788
x0154 =    -1.4315
y0154 =     4.7831
z0154 =    14.3528
x0155 =    -4.6080
y0155 =     4.9660
z0155 =    13.3390
x0156 =    -6.3660
y0156 =     6.5160
z0156 =    13.0800
x0157 =    -5.0378
y0157 =     4.6120
z0157 =    14.2877
x0158 =    -5.1990
y0158 =     6.2880
z0158 =    12.8850
x0159 =    -4.7251
y0159 =     4.1957
z0159 =    12.5625
x0160 =    -3.5370
y0160 =     5.1484
z0160 =    13.5113
x0161 =    -2.1900
y0161 =    22.5140
z0161 =     8.4700
x0162 =    -3.7057
y0162 =    21.0885
z0162 =     7.9189
x0163 =    -2.7240
y0163 =    21.2490
z0163 =     8.0940
x0164 =    -2.7205
y0164 =    22.8206
z0164 =     9.3836
x0165 =    -2.3479
y0165 =    23.2206
z0165 =     7.6419
x0166 =    -1.1091
y0166 =    22.4788
z0166 =     8.6713
x0167 =     4.4810
y0167 =     8.3710
z0167 =    11.6340
x0168 =     6.3600
y0168 =     8.8700
z0168 =    10.1430
x0169 =     4.9075
y0169 =     7.4735
z0169 =    12.1058
x0170 =     5.1340
y0170 =     8.9530
z0170 =    10.3610
x0171 =     4.3710
y0171 =     9.2245
z0171 =    12.3191
x0172 =     3.5262
y0172 =     7.9184
z0172 =    11.3282
x0173 =     2.2990
y0173 =    21.7650
z0173 =     7.8190
x0174 =     4.0730
y0174 =    20.1380
z0174 =     7.7670
x0175 =     2.5881
y0175 =    22.2414
z0175 =     8.7674
x0176 =     2.8600
y0176 =    20.3500
z0176 =     7.7600
x0177 =     2.6846
y0177 =    22.2911
z0177 =     6.9333
x0178 =     1.1994
y0178 =    21.7846
z0178 =     7.7986
 
HrmStr1   * * 331. 1.09
HrmBnd1   * * * 20. 90.
AmbTrs    * * * * 0 0 0 0 0.0 0.0 1.15 0.0 3.0
VDW       MG 1.09 0.25
VDW       NT 1.8240 0.17

   1-0178 0
6-31G*  
****
 
END