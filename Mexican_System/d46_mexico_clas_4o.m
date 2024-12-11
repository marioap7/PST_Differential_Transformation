% Bus data format
% bus: 
% col1 number
% col2 voltage magnitude(pu)
% col3 voltage angle(degree)
% col4 p_gen(pu)
% col5 q_gen(pu),
% col6 p_load(pu)
% col7 q_load(pu)
% col8 G shunt(pu)
% col9 B shunt(pu)
% col10 bus_type
%       bus_type - 1, swing bus
%                - 2, generator bus (PV bus)
%                - 3, load bus (PQ bus)
% col11 q_gen_max(pu)
% col12 q_gen_min(pu)

bus = [ ...
1	1	0	6.000023945	-0.871960072	0	0	0	0	1
2	1	-9.814074074	6.576726772	-1.195682896	0	0	0	0	2
3	1	-8.22173106	5.399999848	-0.782792238	0	0	0	0	2
4	1	-12.2883728	5.99999912	-1.738383145	0	0	0	0	2
5	1	-18.37125205	7.979998869	-0.725773338	0	0	0	0	2
6	1	-21.79480085	3.899998145	-1.318592802	0	0	0	0	2
7	1	-11.57837625	2.399999648	-0.440229732	0	0	0	0	2
8	1	-27.7163444	2.99999918	-0.03222123	0	0	0	0	2
9	1	-27.3973743	1.19999985	-0.330838874	0	0	0	0	2
10	1	-26.64904923	5.399999325	-1.055692338	0	0	0	0	2
11	1	-31.0731453	2.699999284	-0.349728307	0	0	0	0	2
12	1	-26.4081258	0.599999902	-0.296160601	0	0	0	0	2
13	1	-36.65640905	0.899999943	0.425697002	0	0	0	0	2
14	1	-29.93293499	0.599999968	-0.043253152	0	0	0	0	2
15	1	-29.55367005	0.299999984	-0.015076139	0	0	0	0	2
16	1	-26.99878718	5.519999826	-0.396134425	0	0	0	0	2
17	1	-29.93989584	1.199999979	0.031985686	0	0	0	0	2
18	1	-29.43363172	1.139999991	0.014705879	0	0	0	0	2
19	1	-19.28015904	2.999999001	-0.615202266	0	0	0	0	2
20	1	-37.9254479	2.699999899	-0.322014296	0	0	0	0	2
21	1	-34.48454865	0.899999897	-0.03692988	0	0	0	0	2
22	1	-36.36524134	0.599999995	0.023252436	0	0	0	0	2
23	1	-36.94595802	1.439999737	-0.296705697	0	0	0	0	2
24	1	-35.75372191	1.019999899	-0.071119111	0	0	0	0	2
25	1	-35.19412149	3.599999926	0.212481425	0	0	0	0	2
26	1	-3.889929913	1.079999989	-0.180904174	0	0	0	0	2
27	1	-14.20248055	0.839999903	-0.335829138	0	0	0	0	2
28	1	-12.56889261	1.799999813	-0.597740253	0	0	0	0	2
29	1	-6.778047732	0.185999999	-0.125288252	0	0	0	0	2
30	1	-21.96971694	0.545999959	-0.004083767	0	0	0	0	2
31	1	-22.41272895	1.739999769	-0.258939245	0	0	0	0	2
32	1	-29.32419902	1.799999948	0.536279244	0	0	0	0	2
33	1	-5.827126407	1.35	0.461919897	0	0	0	0	2
34	1	-5.827126407	1.35	0.461919897	0	0	0	0	2
35	1	-9.856717505	1.679999128	-1.519088772	0	0	0	0	2
36	1	-26.63813532	5.399999699	-0.986705156	0	0	0	0	2
37	1	-26.06310586	4.199999767	-0.614737323	0	0	0	0	2
38	1	1.758157733	3.479999827	-0.339392006	0	0	0	0	2
39	1	0.146627055	3.479999682	-0.441674794	0	0	0	0	2
40	1	-3.57194918	1.739999981	0.111295329	0	0	0	0	2
41	1	-4.798578899	0.359999996	-0.021183722	0	0	0	0	2
42	1	-3.554188874	0.179999991	0.030346152	0	0	0	0	2
43	1	-23.10339887	0.659999824	-0.242936688	0	0	0	0	2
44	1	-21.94712925	1.199999645	-0.140186842	0	0	0	0	2
45	1	-25.3354446	2.999999939	-0.126371655	0	0	0	0	2
46	1	-26.44597748	1.379998302	-0.253439899	0	0	0	0	2
47	1.043808599	-23.72823131	0	0	0	0	0	0	3
48	1.004933405	-9.976911944	0	0	0.798002074	0.258000088	0	0	3
49	1.008678211	-12.24305019	0	0	0.76805079	0.252015864	0	-1.75	3
50	1.007616378	-13.72146184	0	0	2.286058242	0.750023685	0	0	3
51	1.018509939	-18.36373003	0	0	3.743994245	1.230100981	0	-0.41	3
52	1.02183523	-20.24721868	0	0	8.49E-05	2.56E-05	0	0	3
53	1.043808599	-23.72823131	0	0	-0.000221904	4.60E-05	0	0	3
54	1.056848935	-26.00643639	0	0	-2.92E-05	2.49E-05	0	0	3
55	1.047468802	-28.02993387	0	0	0.000135677	3.18E-05	0	0	3
56	1.055022327	-21.4300524	0	0	0.000158845	4.92E-05	0	0	3
57	1.04304427	-19.45874105	0	0	-1.62E-05	1.83E-05	0	0	3
58	1.039469277	-23.84802614	0	0	1.589938578	0.522010219	0	0	3
59	0.996282884	-32.03036506	0	0	-1.33E-05	9.18E-06	0	-0.82	3
60	1.014202605	-33.61406968	0	0	-1.10E-06	2.29E-06	0	0	3
61	1.008943811	-34.12112859	0	0	0.509998353	0.168000139	0	0	3
62	1.01121722	-33.91785664	0	0	0.629997937	0.20600014	0	0	3
63	1.007384645	-34.23316729	0	0	1.013997375	0.336000166	0	0	3
64	1.036251843	-28.36292973	0	0	4.83E-05	-9.90E-06	0	-0.672	3
65	1.021829333	-30.24074234	0	0	0.533990878	0.17399999	0	0	3
66	1.054934138	-30.43760637	0	0	2.399953803	0.786016011	0	0	3
67	0.964527065	-34.94485127	0	0	1.589998306	0.521999962	0	0	3
68	1.021118498	-30.7650806	0	0	0.59999627	0.198000876	0	0	3
69	1.009172467	-32.5501098	0	0	-1.15E-05	2.83E-05	0	0	3
70	0.991273932	-32.98979098	0	0	1.403991024	0.462002159	0	0	3
71	1.024811246	-28.99584223	0	0	0.972036771	0.329980279	0	-0.672	3
72	1.027078091	-30.62113899	0	0	0.947987079	0.312003995	0	0	3
73	0.996425514	-33.83975473	0	0	-8.74E-06	1.56E-05	0	0	3
74	1.030631258	-27.05050028	0	0	-1.85E-05	3.53E-06	0	0	3
75	1.021337669	-24.65784988	0	0	-3.78E-06	1.04E-05	0	-1.2	3
76	1.020801502	-25.66945564	0	0	-5.57E-07	-3.15E-07	0	0	3
77	1.008077363	-22.09303386	0	0	0.48000881	0.162000824	0	-1.5	3
78	1.010872397	-30.17721213	0	0	-1.22E-05	8.41E-06	0	0	3
79	1.006450097	-31.74718959	0	0	0.485999793	0.162000062	0	0	3
80	1.018764654	-25.79920995	0	0	1.079999095	0.354001292	0	0	3
81	1.01500889	-27.17020473	0	0	0.852000458	0.281999865	0	0	3
82	1.021067215	-27.53248191	0	0	1.37E-07	-7.37E-08	0	0	3
83	1.01371582	-30.05187737	0	0	3.023978073	0.99601117	0	0	3
84	1.023002361	-25.29183647	0	0	0.360011044	0.120002486	0	0	3
85	1.048334856	-27.35513907	0	0	3.12E-06	3.86E-06	0	0	3
86	1.011286837	-22.06874931	0	0	1.139997919	0.372009522	0	-1.888	3
87	1.022807938	-24.27855927	0	0	2.729995733	0.894001244	0	0	3
88	1.015292669	-32.22347993	0	0	-2.57E-06	9.91E-07	0	0	3
89	1.001737378	-30.5315395	0	0	0.180002981	0.059998912	0	-0.75	3
90	1.01173076	-29.64735343	0	0	4.70E-06	7.52E-07	0	0	3
91	1.005191331	-31.07771427	0	0	1.211999453	0.402000165	0	0	3
92	0.999459941	-31.57041316	0	0	-6.30E-07	8.82E-07	0	0	3
93	1.011407889	-33.80720699	0	0	3.371988351	1.110003094	0	0	3
94	1.006498929	-33.98845975	0	0	1.001998062	0.330000103	0	0	3
95	1.011347973	-33.99094917	0	0	0.425998584	0.138000139	0	0	3
96	1.005203913	-33.08962133	0	0	4.445998282	1.619999702	0	0	3
97	1.01540829	-28.11809764	0	0	3.62E-07	-4.57E-07	0	0	3
98	0.996101756	-32.41604585	0	0	1.55E-06	1.27E-07	0	0	3
99	1.001954431	-31.26773061	0	0	6.03E-07	4.16E-08	0	0	3
100	1.006248204	-33.49314145	0	0	2.351997311	0.77400065	0	0	3
101	1.005380715	-33.89829828	0	0	0.47999925	0.156000049	0	0	3
102	0.995523241	-32.97679253	0	0	9.60E-07	-6.13E-07	0	0	3
103	0.992123385	-33.14074069	0	0	3.719989957	1.224000554	0	0	3
104	1.001892989	-34.63544077	0	0	-8.86E-07	3.38E-07	0	0	3
105	0.982613575	-38.89256564	0	0	1.739999259	0.570000007	0	0	3
106	1.007966069	-33.56349168	0	0	1.301997285	0.432001371	0	0	3
107	1.001629626	-29.75389658	0	0	3.10E-06	4.33E-07	0	0	3
108	1.001378655	-26.77746261	0	0	0.228000958	0.078000035	0	0	3
109	1.006831647	-33.95016508	0	0	3.50399426	1.152001966	0	0	3
110	1.006191266	-30.45848195	0	0	2.20E-06	2.92E-07	0	-1.85	3
111	0.999722959	-32.1825004	0	0	3.335999768	0.444000017	0	0	3
112	0.999939326	-30.55722754	0	0	1.241999983	0.359999997	0	0	3
113	1.011907005	-33.07902797	0	0	1.607999338	0.528001308	0	0	3
114	1.017084443	-38.09243428	0	0	2.38E-06	2.48E-06	0	-1.522	3
115	1.023301328	-35.51621163	0	0	1.91E-06	2.19E-06	0	0	3
116	1.013521296	-39.65724474	0	0	1.775997394	0.516000357	0	0	3
117	1.011420458	-39.75323673	0	0	1.775999735	0.516000044	0	0	3
118	1.0194181	-38.38738122	0	0	2.19E-06	3.16E-06	0	-0.772	3
119	1.01312401	-39.65802973	0	0	1.775997482	0.516000232	0	0	3
120	1.012337697	-29.94072965	0	0	1.782005154	0.52200311	0	0	3
121	0.999243878	-37.86263049	0	0	-5.61E-07	3.31E-07	0	-0.7	3
122	1.006398286	-39.3564442	0	0	2.525999117	0.738000179	0	0	3
123	0.994005221	-41.37340437	0	0	1.181999827	0.34200004	0	0	3
124	1.022683041	-41.2629093	0	0	-4.19E-07	1.36E-06	0	0	3
125	1.011323612	-42.47325712	0	0	2.879997345	0.840000373	0	0	3
126	0.997357127	-42.49970594	0	0	1.709999677	0.498000145	0	0	3
127	0.994422281	-38.41081901	0	0	-4.85E-07	-7.09E-08	0	-1.2	3
128	1.003622897	-39.48245519	0	0	1.002000159	0.294000113	0	0	3
129	0.998262713	-38.56403149	0	0	2.790000383	0.93000006	0	0	3
130	0.995861369	-36.61959992	0	0	-8.40E-07	1.35E-07	0	0	3
131	1.004204734	-37.88079073	0	0	3.16199845	0.924000323	0	0	3
132	0.985276676	-41.54286603	0	0	1.505999833	0.438000012	0	0	3
133	0.999250423	-38.52627812	0	0	4.16E-08	-7.23E-10	0	0	3
134	1.003048391	-18.42675533	0	0	1.92E-07	1.16E-08	0	0	3
135	1.001450756	-14.43931788	0	0	0.929999815	0.305999904	0	0	3
136	1.006491604	-5.888446038	0	0	1.23E-07	6.05E-08	0	0	3
137	1.018400281	-16.66069602	0	0	1.338000028	0.444000714	0	0	3
138	1.020335668	-15.74433343	0	0	0.34800072	0.11400005	0	0	3
139	1.024573365	-16.34996203	0	0	0.587999429	0.192000075	0	0	3
140	1.027119524	-15.92086823	0	0	0.353999441	0.113999992	0	0	3
141	1.042600161	-15.11848783	0	0	-9.24E-07	4.34E-07	0	0	3
142	1.013360051	-25.37828583	0	0	1.025999455	0.324000905	0	0	3
143	1.009610115	-25.51485886	0	0	0.233999822	0.07800018	0	0	3
144	1.056868238	-18.21545898	0	0	-1.63E-06	2.04E-06	0	0	3
145	1.002210089	-25.370974	0	0	0.479999988	0.155999993	0	0	3
146	0.988837518	-31.92983159	0	0	-1.34E-07	-4.46E-08	0	0	3
147	1.005455007	-32.51605874	0	0	0.599999814	0.204000088	0	0	3
148	1.000924994	-28.14840518	0	0	0.749999569	0.246000209	0	0	3
149	0.997742903	-28.94490946	0	0	0.299999808	0.102000024	0	0	3
150	0.979642105	-33.24860853	0	0	1.085999882	0.383999996	0	0	3
151	0.984784104	-32.61440024	0	0	0.108000118	0.035999968	0	0	3
152	0.97460646	-33.93439526	0	0	0.995999877	0.33	0	0	3
153	0.99656679	-29.48173495	0	0	0.449999582	0.150000067	0	0	3
154	1.041813809	-14.85764008	0	0	0.581998186	0.192000286	0	0	3
155	1.041998718	-14.77102632	0	0	0.431998704	0.144000045	0	0	3
156	1.03627505	-13.52005596	0	0	0.947996871	0.318000188	0	0	3
157	1.048344651	-12.97214724	0	0	-3.69E-07	1.14E-06	0	0	3
158	1.026194479	-11.8229299	0	0	-3.14E-07	-2.25E-06	0	-1.55	3
159	1.027816028	-11.38124917	0	0	1.50E-06	-3.07E-06	0	-1.2	3
160	1.043301065	-12.81315062	0	0	0.911998146	0.306001354	0	0	3
161	1.009371219	-3.183879807	0	0	2.50E-07	-1.75E-07	0	0	3
162	0.995814177	-5.107867877	0	0	0.833999829	0.275999985	0	0	3
163	1.007977666	-2.821869151	0	0	0.300002802	0.101998182	0	-1.8	3
164	1.009075241	-2.157291288	0	0	0.300000505	0.101999982	0	0	3
165	1.008351094	-2.7967276	0	0	0.408003257	0.192000469	0	-1.544	3
166	1.026807256	-6.657556659	0	0	4.26E-06	-3.60E-06	0	-1.2	3
167	1.016353668	-5.612562787	0	0	0.390000393	0.125999465	0	0	3
168	1.03870726	-10.41831748	0	0	0.88199776	0.287999888	0	0	3
169	1.020663642	-8.570006659	0	0	4.25E-06	-3.81E-06	0	-1.8	3
170	1.041055287	-10.30950909	0	0	0.342000487	0.114001037	0	0	3
171	1.01042626	-5.567914269	0	0	-3.50E-08	-7.24E-08	0	0	3
172	1.007084936	-6.535616918	0	0	0.149999948	0.047999958	0	0	3
173	1.002055443	-7.060777649	0	0	0.545999925	0.180000013	0	0	3
174	1.005165747	-7.212744005	0	0	0.149999986	0.04800003	0	0	3
175	0.99799887	-6.730216321	0	0	0.840000029	0.275999993	0	0	3
176	1.003272667	-7.153269735	0	0	0.228000008	0.078000001	0	0	3
177	0.997941789	-7.212568736	0	0	2.62E-07	-1.73E-07	0	0	3
178	1.008143441	-7.702961379	0	0	0.671999961	0.222000035	0	0	3
179	1.040090219	-12.72610772	0	0	0.659998764	0.222000417	0	0	3
180	0.987789759	-8.020213937	0	0	2.73	0.912	0	0	3
181	0.983135332	-7.930848083	0	0	-1.42E-10	-3.86E-11	0	0	3
182	1.030595423	-18.45078357	0	0	0.395998713	0.132006311	0	-1.7	3
183	1.025711021	-12.39780719	0	0	1.271996141	0.426000238	0	0	3
184	1.042536836	-14.99625065	0	0	0.371998651	0.120002334	0	0	3
185	1.027764599	-12.50288017	0	0	1.53E-06	8.89E-07	0	-0.7	3
186	1.04391389	-12.457867	0	0	0.480002889	0.156000201	0	0	3
187	1.039795103	-16.10436384	0	0	-2.16E-08	7.08E-08	0	0	3
188	1.045561785	-15.27700519	0	0	-1.10E-07	1.19E-08	0	0	3
189	1.021675879	-38.9876036	0	0	-1.31E-07	9.90E-07	0	0	3
190	1.013366509	-40.03717187	0	0	1.77599761	0.516000236	0	0	3];
%
% Line data format
% line: from bus, to bus, resistance(pu), reactance(pu),
%       line charging(pu), tap ratio, phase shift(deg)
line = [ ...
48	49	0.0007	0.0086	0.9714	0	0
49	50	0.0005	0.0063	0.7192	0	0
50	52	0.0028	0.0338	0.9878	0	0
52	51	0	-0.0101	0	0	0
50	52	0.0028	0.0338	0.9878	0	0
52	51	0	-0.0101	0	0	0
50	52	0.0028	0.0338	0.9878	0	0
52	51	0	-0.0101	0	0	0
49	56	0.0022	0.0266	3.3358	0	0
57	53	0.0014	0.0177	2.1908	0	0
57	56	0	-0.0083	0	0	0
51	55	0.0041	0.0552	1.444	0	0
53	55	0	-0.0259	0	0	0
53	69	0.0042	0.0574	1.4215	0	0
51	55	0.0041	0.0552	1.444	0	0
53	55	0	-0.0259	0	0	0
53	69	0.0042	0.0574	1.4215	0	0
70	83	0.0062	0.0381	0.2875	0	0
69	64	0	-0.027	0	0	0
69	64	0	-0.027	0	0	0
53	73	0.0038	0.0531	1.3141	0	0
71	73	0	-0.0249	0	0	0
71	74	0	-0.0124	0	0	0
71	64	0.0008	0.092	0.2688	0	0
74	59	0.0025	0.0317	0.9003	0	0
59	92	0.0003	0.0039	0.457	0	0
64	78	0.0017	0.023	0.591	0	0
64	78	0.0017	0.023	0.591	0	0
75	89	0.0041	0.053	1.5312	0	0
75	86	0.002	0.0263	3.035	0	0
75	84	0.0028	0.0345	1.008	0	0
75	77	0.0014	0.0176	0.5009	0	0
77	78	0.0048	0.0577	1.596	0	0
77	78	0.0048	0.0577	1.596	0	0
89	78	0.0015	0.0181	0.5429	0	0
89	78	0.0015	0.0181	0.5429	0	0
87	80	0.0315	0.1113	0.4645	0	0
80	76	0.0007	0.0046	0.0088	0	0
80	81	0.0097	0.0642	0.121	0	0
81	68	0.0168	0.1116	0.2105	0	0
68	65	0.007	0.0465	0.0967	0	0
68	72	0.0176	0.1179	0.2233	0	0
68	88	0.0126	0.0837	0.1579	0	0
68	96	0.0151	0.1004	0.1894	0	0
97	96	0.0072	0.0567	0.4608	0	0
97	76	0.0068	0.0455	0.3446	0	0
82	76	0.0124	0.0836	0.1582	0	0
85	83	0.0056	0.0372	0.2806	0	0
83	54	0.0081	0.0497	0.3751	0	0
82	81	0.0049	0.0326	0.0614	0	0
82	68	0.0154	0.1024	0.193	0	0
84	64	0.0042	0.0513	1.4984	0.00E+00	0.00E+00
84	71	0.004	0.0508	1.4369	0	0.00E+00
182	86	0.0016	0.0234	2.4142	0.00E+00	0.00E+00
182	185	0.0046	0.0584	1.6584	0	0.00E+00
185	184	0.002	0.0253	0.7175	0	0.00E+00
184	182	0.0037	0.0465	1.3199	0.00E+00	0.00E+00
78	92	0.0004	0.0044	0.5286	0	0
89	98	0.0005	0.0062	0.7374	0.00E+00	0.00E+00
90	109	0.0105	0.0698	0.1315	0.00E+00	0.00E+00
90	131	0.0213	0.1263	0.2354	0	0
88	96	0.0015	0.0088	0.0735	0	0
88	79	0.0046	0.0347	0.0702	0	0
91	79	0.0022	0.0165	0.0334	0.00E+00	0.00E+00
91	90	0.0021	0.0135	0.1052	0	0
93	96	0.0043	0.0344	0.0691	0	0
93	62	0.0015	0.012	0.027	0	0
93	95	0.0037	0.0305	0.0735	0	0
93	94	0.0013	0.0104	0.023	0.00E+00	0.00E+00
94	96	0.003	0.0242	0.054	0	0
96	131	0.0308	0.1771	0.3586	0	0
98	102	0.0002	0.0027	0.32	0	0
100	96	0.0013	0.0102	0.092	0.00E+00	0.00E+00
100	109	0.0025	0.0194	0.04	0.00E+00	0.00E+00
100	101	0.002	0.013	0.0246	0.00E+00	0.00E+00
101	109	0.0015	0.0102	0.0193	0.00E+00	0.00E+00
67	107	0.0085	0.0567	0.428	0	0
107	106	0.0085	0.0567	0.428	0.00E+00	0.00E+00
107	108	0.0015	0.0187	0.2327	0	0
106	72	0.0181	0.1114	0.2103	0	0
106	60	0.0037	0.0251	0.1899	0	0
62	60	0.0016	0.0123	0.0257	0.00E+00	0.00E+00
61	60	0.0015	0.0114	0.023	0	0
95	60	0.0029	0.0219	0.045	0	0
65	72	0.0023	0.0139	0.1052	0.00E+00	0.00E+00
106	104	0.0092	0.0599	0.1131	0	0
104	109	0.0036	0.0242	0.0456	0	0
63	109	0.0019	0.0148	0.0298	0.00E+00	0.00E+00
63	61	0.0008	0.0065	0.013	0	0
63	60	0.0043	0.033	0.0674	0.00E+00	0.00E+00
110	113	0.0023	0.0292	3.3168	0	0
110	111	0.0011	0.0143	0.4061	0.00E+00	0.00E+00
113	102	0.0009	0.0113	1.286	0	0
103	102	0.0004	0.0042	0.128	0	0
59	103	0.0006	0.0065	0.1993	0	0
110	115	0.0038	0.0524	1.4892	0	0
114	115	0.0017	0.0215	0.6092	0.00E+00	0.00E+00
118	114	0.0005	0.0062	0.1754	0.00E+00	0.00E+00
118	124	0.0039	0.0501	1.4167	0.00E+00	0.00E+00
118	189	0.0013	0.0169	0.4738	0	0
118	120	0.0039	0.0503	1.4166	0	0
114	120	0.0044	0.0569	1.5841	0.00E+00	0.00E+00
120	115	0.0034	0.0441	1.2455	0	0
121	115	0.002	0.025	0.7107	0.00E+00	0.00E+00
127	121	0.0023	0.0298	0.8461	0	0
116	117	0.0015	0.0097	0.313	0	0
116	190	0.0017	0.0116	0.3516	0.00E+00	0.00E+00
114	189	0.0009	0.0119	0.3384	0	0
189	124	0.0038	0.0429	1.2184	0	0
118	127	0.0043	0.0549	1.5516	0.00E+00	0.00E+00
119	117	0.0009	0.0058	0.1758	0	0
127	122	0.0173	0.1161	0.2198	0	0
119	122	0.0246	0.1653	0.313	0	0
127	130	0.0015	0.0208	0.5869	0.00E+00	0.00E+00
130	89	0.0031	0.0394	1.1131	0.00E+00	0.00E+00
128	126	0.0056	0.0372	0.2806	0	0
128	123	0.0035	0.0233	0.0438	0	0
123	126	0.0109	0.0725	0.1368	0.00E+00	0.00E+00
126	125	0.0076	0.0502	0.3788	0	0
125	129	0.0104	0.0688	0.5192	0	0
128	131	0.0134	0.0785	0.1427	0.00E+00	0.00E+00
128	132	0.007	0.0411	0.0747	0	0
131	129	0.013	0.0865	0.6525	0	0
133	132	0.0132	0.0882	0.167	0.00E+00	0.00E+00
117	190	0.0035	0.0232	0.1758	0	0
136	134	0.0298	0.2023	0.3888	0	0
134	135	0	-0.0668	0	0.00E+00	0.00E+00
135	138	0.0257	0.1735	0.3369	0	0
138	137	0.0018	0.0121	0.0914	0	0
158	141	0.0059	0.0756	2.1458	0.00E+00	0.00E+00
138	140	0.0016	0.011	0.0851	0	0
139	140	0.0019	0.0128	0.0248	0	0
157	155	0.0032	0.0529	0.1002	0.00E+00	0.00E+00
137	139	0.0035	0.0233	0.0408	0.00E+00	0.00E+00
139	188	0.0146	0.0973	0.1834	0	0
188	187	0	-0.0668	0	0.00E+00	0.00E+00
137	142	0.0204	0.1358	1.0243	0	0
143	148	0.0055	0.0371	0.2813	0	0
143	148	0.0055	0.0371	0.2813	0	0
143	144	0.0241	0.1606	0.3049	0	0
148	153	0.004	0.0269	0.051	0.00E+00	0.00E+00
142	143	0.0047	0.0316	0.2393	0	0
146	153	0.0063	0.0423	0.32	0	0
148	149	0.0024	0.0158	0.0299	0.00E+00	0.00E+00
149	153	0.0024	0.0158	0.0299	0	0
146	151	0.005	0.0302	0.2307	0.00E+00	0.00E+00
151	150	0.001	0.0068	0.0512	0	0
151	152	0.0079	0.0523	0.1042	0	0
152	150	0.0032	0.0214	0.0404	0	0
157	154	0.009	0.0604	0.1143	0	0
154	155	0.0011	0.0074	0.0141	0	0
154	187	0.0146	0.0993	0.1925	0	0
157	160	0.0017	0.0112	0.0842	0	0
158	159	0.0005	0.0057	0.1624	0	0
185	183	0.0007	0.0093	0.264	0	0
158	183	0.0006	0.0072	0.2031	0	0
169	158	0.0017	0.0207	2.4192	0.00E+00	0.00E+00
170	144	0.0132	0.1687	0.5237	0.00E+00	0.00E+00
170	168	0.001	0.007	0.0527	0.00E+00	0.00E+00
163	169	0.0019	0.0229	2.6746	0	0
164	161	0.0021	0.0189	0.1055	0.00E+00	0.00E+00
161	167	0.0052	0.0344	0.2596	0	0
167	168	0.0155	0.1032	0.1947	0	0
168	160	0.0259	0.1721	0.3245	0	0
164	171	0.0119	0.0791	0.5964	0	0
164	161	0.0029	0.0195	0.0369	0.00E+00	0.00E+00
172	173	0.0049	0.0179	0.0044	0	0
172	174	0.0785	0.2867	0.0712	0	0
173	174	0.0761	0.3091	0.076	0.00E+00	0.00E+00
181	174	0.013	0.434	0.1116	0	0
178	176	0.005	0.0295	0.0318	0.00E+00	0.00E+00
177	175	0.0017	0.0111	0.0844	0	0
174	178	0.0366	0.142	0.152	0	0
186	177	0.0319	0.212	0.3999	0	0
185	159	0.0006	0.0072	0.2031	0	0
160	186	0.0034	0.0065	0.0144	0	0
165	166	0.0011	0.0139	3.5526	0.00E+00	0.00E+00
165	163	0.0001	0.0007	0.0203	0	0
159	166	0.0014	0.0178	2.0304	0	0
185	159	0.003	0.0384	1.0898	0	0
179	186	0.0045	0.0085	0.194	0.00E+00	0.00E+00
53	47	0	0.001	0	0	0
48	3	0	0.0057	0	0	0
49	2	0	0.0065	0	0	0
50	4	0	0.0042	0	0.00E+00	0.00E+00
50	7	0	0.0157	0	0	0
56	58	0	0.0291	0	0.00E+00	0.00E+00
53	54	0	0.0275	0	9.75E-01	0.00E+00
59	60	0	0.0145	0	9.70E-01	0.00E+00
64	65	0	0.0291	0	1.025	0
64	66	0	0.017	0	0.97	0
75	76	0	0.0095	0	0	0
77	5	0	0.0082	0	0	0
84	85	0	0.0275	0	0.96	0
86	87	0	0.0137	0	0.98	0
80	43	0	0.0726	0	0	0
81	44	0	0.077	0	0	0
83	46	0	0.0462	0	0	0
84	6	0	0.016	0	0	0
86	19	0	0.0164	0	0	0
78	88	0	0.0275	0	0.97	0
89	90	0	0.0065	0	0.98	0
89	8	0	0.0164	0	0	0
90	9	0	0.0331	0	0	0
90	10	0	0.0098	0	0	0
92	93	0	0.0111	0	0.97	0
96	11	0	0.0131	0	0	0
97	12	0	0.0505	0	0	0
105	13	0	0.0426	0	0	0
104	105	0	0.087	0	0	0
98	100	0	0.0111	0	0.98	0
100	99	0	0.0435	0	0	0
99	14	0	0.0389	0	0	0
99	15	0	0.0999	0	0	0
102	109	0	0.0062	0	0.98	0
108	45	0	0.0084	0	0	0
71	72	0	0.0291	0	0.98	0
110	16	0	0.011	0	0	0
111	17	0	0.0326	0	0	0
110	112	0	0.017	0	0	0
112	18	0	0.0172	0	0	0
114	116	0	0.011	0	0	0
121	122	0	0.019	0	0.98	0
122	24	0	0.062	0	0	0
118	119	0	0.011	0	0	0
189	190	0	0.019	0	0	0
127	128	0	0.011	0	0.98	0
130	131	0	0.019	0	0.98	0
128	20	0	0.0101	0	0	0
133	22	0	0.0628	0	0	0
117	23	0	0.0344	0	0	0
129	25	0	0.0163	0	0	0
124	125	0	0.011	0	0	0
131	21	0	0.0661	0	0	0
136	26	0	0.0325	0	0	0
141	140	0	0.019	0	0	0
138	28	0	0.0314	0	0	0
137	27	0	0.052	0	0	0
151	32	0	0.0314	0	0	0
146	147	0	0.0173	0	0.98	0
143	145	0	0.0385	0	0	0
145	30	0	0.1089	0	0	0
143	31	0	0.0314	0	0	0
158	156	0	0.0339	0	0.98	0
169	170	0	0.021	0	0.98	0
158	157	0	0.021	0	0.97	0
181	180	0	0.051	0	0.99	0
171	172	0	0.0231	0	0	0
29	174	0	0.041	0	0	0
177	178	0	0.0268	0	0.98	0
176	41	0	0.1145	0	0	0
162	42	0	0.15	0	0	0
161	162	0	0.0516	0	0	0
175	40	0	0.0316	0	0	0
175	176	0	0.046	0	0	0
164	38	0	0.0198	0	0	0
163	39	0	0.015	0	0	0
165	1	0	0.0082	0	0	0
163	164	0	0.021	0	0	0
180	33	0	0.028	0	0	0
180	34	0	0.028	0	0	0
185	186	0	0.0104	0	0.97	0
186	35	0	0.0282	0	0	0
120	36	0	0.0108	0	0	0
120	37	0	0.0163	0	0	0];

% Machine data format
%       1. machine number,
%       2. bus number,
%       3. base mva,
%       4. leakage reactance x_l(pu),
%       5. resistance r_a(pu),
%       6. d-axis sychronous reactance x_d(pu),
%       7. d-axis transient reactance x'_d(pu),
%       8. d-axis subtransient reactance x"_d(pu),
%       9. d-axis open-circuit time constant T'_do(sec),
%      10. d-axis open-circuit subtransient time constant
%                T"_do(sec),
%      11. q-axis sychronous reactance x_q(pu),
%      12. q-axis transient reactance x'_q(pu),
%      13. q-axis subtransient reactance x"_q(pu),
%      14. q-axis open-circuit time constant T'_qo(sec),
%      15. q-axis open circuit subtransient time constant
%                T"_qo(sec),
%      16. inertia constant H(sec),
%      17. damping coefficient d_o(pu),
%      18. dampling coefficient d_1(pu),
%      19. bus number
%

mac_con = ...
  [1 1 100 0 0 0.094373401534526849 0.013938618925831203 0 5.2 0 0.093350383631713552 ...
   0.024296675191815855 0 0.47 0 53.958 18.85 0 1 0 0 1 1;
   2 2 100 0 0 0.045612923661704149 0.015837820715869498 0.036 7.4 0.029 ...
   0.025340513145391194 0.025340513145391194 0.036 1000 0.034 73.08455 ...
   18.85 0 2 0 0 1 1;
   3 3 100 0 0 0.078534031413612565 0.025130890052356018 0.034 5.2 0.029 ...
   0.0450261780104712 0.0450261780104712 0.034 1000 0.034 41.065 18.85 ...
   0 3 0 0 1 1;
   4 4 100 0 0 0.062385321100917435 0.01908256880733945 0.03 5.53 0.084 ...
   0.034862385321100912 0.034862385321100912 0.03 1000 0.198 52.974000000000004 ...
   18.85 0 4 0 0 1 1;
   5 5 100 0 0 0.094373401534526849 0.013938618925831203 0.0116 5.2 0.032 ...
   0.093350383631713552 0.024296675191815855 0.0243 0.47 0.06 53.958 18.85 ...
   0 5 0 0 1 1;
   6 6 100 0 0 0.21066666666666667 0.042800000000000005 0.262 6.5 0.06 ...
   0.19866666666666666 0.11266666666666666 0.0349 0.7 0.06 44.025 18.85 ...
   0 6 0 0 1 1;
   7 7 100 0 0 0.196 0.05 0.05 6 0.06 0.108 0.108 0.05 1000 0.09 16.15 ...
   18.85 0 7 0 0 1 1;
   8 8 100 0 0 0.239564961787184 0.038212815990593771 0.11 5.5 0.02 0.2278071722516167 ...
   0.066137566137566134 0.0323 0.6 0.03 21.092399999999998 18.85 0 8 0 ...
   0 1 1;
   9 9 100 0 0 0.09 0.09 0 2 0 0.09 0.09 0 0 0 12.6 18.85 0 9 0 0 1 1;
   10 10 100 0 0 0.1596119929453263 0.025455614344503236 0.0733 5.5 0.02 ...
   0.15167548500881836 0.04409171075837743 0 0.6 0.03 31.6386 18.85 0 10 ...
   0 0 1 1;
   11 11 100 0 0 0.31150341685649208 0.035876993166287015 0.0503 5.7 0.06 ...
   0.28246013667425968 0.13268792710706151 0.0286 0.7 0.06 16.64688 18.85 ...
   0 11 0 0 1 1;
   12 12 100 0 0 0.183 0.183 0 2 0 0.183 0.183 0 0 0 6.3 18.85 0 12 0 0 ...
   1 1;
   13 13 100 0 0 0.1947 0.1947 0 2 0 0.1947 0.1947 0 0 0 4.8 18.85 0 13 ...
   0 0 1 1;
   14 14 100 0 0 0.08 0.08 0 2 0 0.08 0.08 0 0 0 7.5 18.85 0 14 0 0 1 1 ...
   ;
   15 15 100 0 0 0.33140000000000003 0.33140000000000003 0 2 0 0.33140000000000003 ...
   0.33140000000000003 0 0 0 3.2 18.85 0 15 0 0 1 1;
   16 16 100 0 0 0.072222222222222215 0.022777777777777782 0.0345 5.88 ...
   0.06 0.045 0.045 0.0345 1000 0.06 36.180000000000007 18.85 0 16 0 0 ...
   1 1;
   17 17 100 0 0 0.41876606683804629 0.057840616966580979 0.181 5.06 0.019 ...
   0.40102827763496146 0.16015424164524422 0 0.6 0.069 12.603600000000002 ...
   18.85 0 17 0 0 1 1;
   18 18 100 0 0 0.25 0.087500000000000008 0.055 5.5 0.06 0.15 0.15 0.055 ...
   1000 0.06 11.744000000000002 18.85 0 18 0 0 1 1;
   19 19 100 0 0 0.22819767441860464 0.040697674418604654 0 5.5 0 0.22529069767441862 ...
   0.22529069767441862 0 1000 0 21.7408 18.85 0 19 0 0 1 1;
   20 20 100 0 0 0.11899999999999998 0.0191 0.0161 5.5 0.02 0.114 0.033 ...
   0.0161 0.6 0.03 42.18 18.85 0 20 0 0 1 1;
   21 21 100 0 0 0.17059999999999997 0.17059999999999997 0 2 0 0.17059999999999997 ...
   0.17059999999999997 0 0 0 4.14 18.85 0 21 0 0 1 1;
   22 22 100 0 0 0.028899999999999995 0.028899999999999995 0 2 0 0.028899999999999995 ...
   0.028899999999999995 0 0 0 4.17 18.85 0 22 0 0 1 1;
   23 23 100 0 0 0.066 0.066 0 2 0 0.066 0.066 0 0 0 11.06 18.85 0 23 0 ...
   0 1 1;
   24 24 100 0 0 0.09 0.09 0 2 0 0.09 0.09 0 0 0 12.6 18.85 0 24 0 0 1 ...
   1;
   25 25 100 0 0 0.20971867007672634 0.033887468030690537 0.1135 5.5 0.032 ...
   0.12915601023017903 0.05754475703324808 0.1135 0.47 0.06 28.6994 18.85 ...
   0 25 0 0 1 1;
   26 26 100 0 0 0.479129923574368 0.076425631981187542 0.22 5.5 0.02 0.4556143445032334 ...
   0.13227513227513227 0 0.6 0.03 10.546199999999999 18.85 0 26 0 0 1 1 ...
   ;
   27 27 100 0 0 0.94574468085106389 0.17872340425531916 0.087 6.5 0.035 ...
   0.92340425531914894 0.37234042553191488 0.087 0.96 0.06 12.708799999999998 ...
   18.85 0 27 0 0 1 1;
   28 28 100 0 0 0.53205849268841388 0.06608548931383576 0 5.75 0 0.43307086614173224 ...
   0.12373453318335208 0 0.46 0 11.236960000000002 18.85 0 28 0 0 1 1;
   29 29 100 0 0 2.4 2.4 0 2 0 2.4 2.4 0 0 0 1.17 18.85 0 29 0 0 1 1;
   30 30 100 0 0 0.1704 0.1704 0 2 0 0.1704 0.1704 0 0 0 6.21 18.85 0 30 ...
   0 0 1 1;
   31 31 100 0 0 0.9 0.0825 0 4.8 0 0.875 0.235 0 0.5 0 5.6 18.85 0 31 ...
   0 0 1 1;
   32 32 100 0 0 0.51252847380410027 0.036731207289293855 0.054 4.8 0.035 ...
   0.49829157175398636 0.13382687927107062 0 0.5 0.07 12.36224 18.85 0 ...
   32 0 0 1 1;
   33 33 100 0 0 0.60227272727272729 0.065454545454545446 0.045 6 0.05 ...
   0.56818181818181823 0.15113636363636365 0.045 0.5 0.06 10.56 18.85 0 ...
   33 0 0 1 1;
   34 34 100 0 0 0.59578947368421042 0.10526315789473684 0.076 3.75 0.03 ...
   0.56105263157894736 0.21052631578947367 0.076 0.4 0.04 14.9055 18.85 ...
   0 34 0 0 1 1;
   35 35 100 0 0 0.49746707193515705 0.061803444782168183 0.028 5.25 0.023 ...
   0.47011144883485312 0.10435663627152987 0.028 0.442 0.055 25.97784 18.85 ...
   0 35 0 0 1 1;
   36 36 100 0 0 0.1596119929453263 0.025455614344503236 0.0733 5.5 0.02 ...
   0.15167548500881836 0.04409171075837743 0.0215 0.6 0.03 31.6386 18.85 ...
   0 36 0 0 1 1;
   37 37 100 0 0 0.2120822622107969 0.027634961439588688 0 4 0 0.20437017994858611 ...
   0.052699228791773779 0 0.52 0 25.207200000000004 18.85 0 37 0 0 1 1;
   38 38 100 0 0 0.18371546149323928 0.030276308054085831 0.092 6.2 0.05 ...
   0.1616696061140506 0.058788947677836566 0.092 0.5 0.06 21.092399999999998 ...
   18.85 0 38 0 0 1 1;
   39 39 100 0 0 0.18371546149323928 0.030276308054085831 0.092 6.2 0.05 ...
   0.1616696061140506 0.058788947677836566 0.092 0.5 0.06 21.092399999999998 ...
   18.85 0 39 0 0 1 1;
   40 40 100 0 0 0.46531791907514453 0.076300578034682084 0.227 5.5 0.05 ...
   0.45664739884393063 0.11560693641618497 0.227 1 0.05 12.213799999999999 ...
   18.85 0 40 0 0 1 1;
   41 41 100 0 0 0.28977272727272729 0.28977272727272729 0 2 0 0.28977272727272729 ...
   0.28977272727272729 0 0 0 3.52 18.85 0 41 0 0 1 1;
   42 42 100 0 0 2.5 0.73529411764705888 0.22 4.7 0.022 1.911764705882353 ...
   1.911764705882353 0.22 1000 0.04 1.0131999999999999 18.85 0 42 0 0 1 ...
   1;
   43 43 100 0 0 0.18 0.18 0 2 0 0.18 0.18 0 0 0 6.3 18.85 0 43 0 0 1 1 ...
   ;
   44 44 100 0 0 0.21034482758620687 0.074137931034482754 0.0307 5 0.05 ...
   0.1586206896551724 0.1586206896551724 0.0307 1000 0.05 9.395999999999999 ...
   18.85 0 44 0 0 1 1;
   45 45 100 0 0 0.13090909090909089 0.034363636363636367 0.063 7.08 0.1 ...
   0.085909090909090907 0.085909090909090907 0.063 1000 0.15 28.71 18.85 ...
   0 45 0 0 1 1;
   46 46 100 0 0 0.62817796610169485 0.10444915254237287 0.058 5.79 0.035 ...
   0.61228813559322026 0.24682203389830507 0.058 0.96 0.06 19.144320000000004 ...
   18.85 0 46 0 0 1 1];
mac_con(:,8) = mac_con(:,13);
% %


% % %
exc_con = [0  1 0  150.000    0.030 0 0   20.000   -20.000
0  2 0  100.000    0.050 0 0   20.000   -20.000
0  3 0  100.000    0.050 0 0   20.000   -20.000
0  4 0  100.000    0.050 0 0   20.000   -20.000
0  5 0  100.000    0.050 0 0   20.000   -20.000
0  6 0  175.000    0.030 0 0   20.000   -20.000
0  7 0  100.000    0.050 0 0   20.000   -20.000
0  8 0  100.000    0.050 0 0   20.000   20.000
0  9 0  100.000    0.050 0 0   20.000   -20.000
0 10 0  100.000    0.050 0 0   20.000   -20.000
0 11 0  100.000    0.050 0 0   20.000   20.000
0 12 0  100.000    0.050 0 0   20.000   20.000
0 13 0  100.000    0.050 0 0   20.000   20.000
0 14 0  100.000    0.050 0 0   20.000   20.000
0 15 0  100.000    0.050 0 0   20.000   -20.000
0 16 0   75.000    0.040 0 0   20.000   -20.000
0 17 0  100.000    0.050 0 0   20.000   -20.000
0 18 0  150.000    0.040 0 0   20.000   20.000
0 19 0  100.000    0.050 0 0   20.000   -20.000
0 20 0  100.000    0.050 0 0   20.000   20.000
0 21 0  100.000    0.050 0 0   20.000   20.000
0 22 0  100.000    0.050 0 0   20.000   20.000
0 23 0  100.000    0.050 0 0   20.000   20.000
0 24 0  100.000    0.050 0 0   20.000   -20.000
0  25 0  100.000    0.050 0 0   20.000  -20.000
0 26 0  100.000    0.050 0 0   20.000   -20.000
0 27 0  100.000    0.050 0 0   20.000   -20.000
0 28 0  100.000    0.050 0 0   20.000   20.000
0 29 0  100.000    0.050 0 0   20.000   20.000
0 30 0  100.000    0.050 0 0   20.000   -20.000
0 31 0  100.000    0.050 0 0   20.000   -20.000
0 32 0  100.000    0.050 0 0   20.000   -20.000
0 33 0  100.000    0.050 0 0   20.000   -20.000
0 34 0  100.000    0.050 0 0   20.000   -20.000
0 35 0  100.000    0.050 0 0   20.000   -20.000
0 36 0  100.000    0.050 0 0   20.000   -20.000
0 37 0  100.000    0.050 0 0   20.000   -20.000
0 38 0   75.000    0.040 0 0   20.000   -20.000
0 39 0  100.000    0.050 0 0   20.000   -20.000
0 40 0  100.000    0.050 0 0   20.000   -20.000
0 41 0  100.000    0.050 0 0   20.000   20.000
0 42 0  100.000    0.050 0 0   20.000   20.000
0 43 0  100.000    0.050 0 0   20.000   20.000
0 44 0  100.000    0.050 0 0   20.000   -20.000
0 45 0  100.000    0.040 0 0   20.000   20.000
0 46 0  100.000    0.050 0 0   20.000   -20.000];
exc_con(:,4)=75; exc_con(:,5)=0.015;
% exc_con(:,9) = -20.00;
% exc_con = [];
%Switching file defines the simulation control
% row 1 col1  simulation start time (s) (cols 2 to 6 zeros)
%       col7  initial time step (s)
% row 2 col1  fault application time (s)
%       col2  bus number at which fault is applied
%       col3  bus number defining far end of faulted line
%       col4  zero sequence impedance in pu on system base
%       col5  negative sequence impedance in pu on system base
%       col6  type of fault  - 0 three phase
%                            - 1 line to ground
%                            - 2 line-to-line to ground
%                            - 3 line-to-line
%                            - 4 loss of line with no fault
%                            - 5 loss of load at bus
%                            - 6 no action
%       col7  time step for fault period (s)
% row 3 col1  near end fault clearing time (s) (cols 2 to 6 zeros)
%       col7  time step for second part of fault (s)
% row 4 col1  far end fault clearing time (s) (cols 2 to 6 zeros)
%       col7  time step for fault cleared simulation (s)
% row 5 col1  time to change step length (s)
%       col7  time step (s)
%
%
%
% row n col1 finishing time (s)  (n indicates that intermediate rows may be inserted)
tao = 1/60;
sw_con = [...
0     0    0    0    0    0    tao;%sets intitial time step
5*tao   144    143  0    0    0    tao; %3 ph fault fault at bus 3
10*tao  0    0    0    0    0    tao; %
% 0.41 0    0    0    0    0    0.01; %clear remote end
6.0  0    0    0    0    0    0]; % end simulation

% ibus_con = [1 0 0 0];

% mac_con(43,:) = [];
% mac_con(41,:) = [];
% mac_con(37,:) = [];
% mac_con(32,:) = [];
% mac_con(31,:) = [];
% mac_con(30,:) = [];
% mac_con(29,:) = [];
% mac_con(28,:) = [];
% mac_con(26,:) = [];
% mac_con(24,:) = [];
% mac_con(23,:) = [];
% mac_con(22,:) = [];
% mac_con(21,:) = [];
% mac_con(19,:) = [];
% mac_con(17,:) = [];
% mac_con(15,:) = [];
% mac_con(14,:) = [];
% mac_con(13,:) = [];
% mac_con(12,:) = [];
% mac_con(10,:) = [];
% mac_con(9,:) = [];
% mac_con(1,:) = [];


