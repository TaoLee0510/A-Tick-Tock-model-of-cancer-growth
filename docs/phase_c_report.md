# Phase C report

Phase C aligns ABM and PDE with option (a), a shared VEGF-guided lattice-tip law,
under `shared_vegf_lattice_v2` and structured schema 15. Published vascular and
tumour models retain their operators, fingerprints, output columns and checkpoints.

The shared law specifies VEGF production/diffusion/decay, guided hops, Yule
branching, local connection termination, hypoxic Poisson roots, unique centerline
growth and vessel-union perfusion. The ABM samples the law; the PDE advances its
conditional first moments with an explicit independent-arrival edge closure.
There is one owned vascular process per standalone model. New overlapping
vessels delete cells and record the six phenotype/stage mass counters.

Parallel work writes disjoint voxels, reductions and individual commits use a
canonical order, work buffers are reused, and VEGF/tip/source support determines
the active range. Caller-supplied consumer bounds avoid full-grid source scans.
The checkpoint includes edge occupancies, tip UIDs/RNG counters, fields and diagnostics.

## Tests and exact arithmetic

| Build | CTest passed | Validation tests | Full-suite seconds |
| --- | ---: | ---: | ---: |
| release | 56/56 | 11 | 2798.3099999999999 |
| hdf5 | 59/59 | 11 | 2798.1199999999999 |

New assert tests cover 2D/3D gradient flux, branching/connection moments, tip
mass accounting, unique edges, active bounds, stale-source clearing, thread
determinism, and restart. The coupling test checks individual deletion and
footprint demand with native and HDF5 checkpoints. The statistical unit test
rejects late undefined zero-reference tests and wrong geometry/time/length models.
All new C++ tests are compiled with `-UNDEBUG`.

Every 416 published native JSON records and every 96 schema-14 trace records
are exactly unchanged in both builds. All 48 new vascular native reports are
exactly equal across the two builds and to the pre-bounds-optimization reports,
including field checksums. There are no build/link warnings in the checked logs.

The existing validation tables are numerically unchanged and retained at full
precision in [the phase B numeric report](phase_b_numeric_results.json).
New tables, all native vascular realizations and raw benchmark samples are in
[the phase C numeric report](phase_c_numeric_results.json).

## Prespecified controlled ensemble

The grid is 256-square, duration 240 hours, ABM A/B groups have 16 disjoint
seeds each, and samples are taken every 24 hours. Tumour migration, activation,
division and exchange are disabled to isolate the vascular mechanism. Nutrient
consumption and vessel perfusion remain coupled. Vessel exclusion is disabled
in this experiment so a fixed consumer realization is preserved.

All post-initialization paired TOST tests pass at the margins recorded before
running: 0.75 relative for centerline, volume and counts, 0.15 absolute for
lesion perfusion. Initialization structural zeros are exact setup checks, not
inferred equivalence. These are broad exploratory mechanism-screening margins.
The result does not demonstrate accurate vessel morphology or biological
interchangeability. The systematic PDE overestimation is visible below.

| Hours | Metric | ABM mean | PDE mean | ABM baseline difference | Model/baseline RMS | Error | Margin | Paired t | Paired p | TOST lower p | TOST upper p | Result |
| ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| 0 | vascular_length | 0 | 0 | 0 | null | null | 0.75 | null | 1 | null | null | INITIALIZATION |
| 0 | perfused_volume | 0 | 0 | 0 | null | null | 0.75 | null | 1 | null | null | INITIALIZATION |
| 0 | lesion_perfused_fraction | 0 | 0 | 0 | null | 0 | 0.14999999999999999 | null | 1 | 0 | 0 | PASS |
| 0 | vascular_branches | 0 | 0 | 0 | null | null | 0.75 | null | 1 | null | null | INITIALIZATION |
| 0 | vascular_anastomoses | 0 | 0 | 0 | null | null | 0.75 | null | 1 | null | null | INITIALIZATION |
| 0 | vascular_roots | 0 | 0 | 0 | null | null | 0.75 | null | 1 | null | null | INITIALIZATION |
| 24 | vascular_length | 409.9375 | 457.80460240610341 | -8.6875 | 1.4891812759193155 | 0.11676683008044739 | 0.75 | 10.401686751481668 | 2.9689841086513843e-08 | 1.2602214944004283e-29 | 1.5615410597187239e-15 | PASS |
| 24 | perfused_volume | 206.15024376529516 | 234.08097356018101 | -3.8801179148965321 | 2.0624282858825289 | 0.1354872508745871 | 0.75 | 14.021740660544847 | 5.0103929632110569e-10 | 1.479385695971675e-30 | 2.4867492481993171e-16 | PASS |
| 24 | lesion_perfused_fraction | 0.3438789498763038 | 0.39046967217834699 | -0.0060508164343887443 | 2.0455819951300342 | 0.046590722302043169 | 0.14999999999999999 | 14.093859971140578 | 4.6630658425535206e-10 | 1.5810225695059485e-19 | 2.2404108779947431e-15 | PASS |
| 24 | vascular_branches | 362.75 | 441.39900233878188 | -1.5625 | 1.0538948884624462 | 0.21681323870098376 | 0.75 | 7.6284482003992364 | 1.5389866522329996e-06 | 1.7206102010467691e-24 | 9.422624868057642e-09 | PASS |
| 24 | vascular_anastomoses | 323.5625 | 325.91587757925862 | -4.875 | 0.03499857354972942 | 0.0072733322905424594 | 0.75 | 0.23292023044017773 | 0.81897176750812739 | 4.304308759032426e-22 | 4.1095375195129292e-10 | PASS |
| 24 | vascular_roots | 114 | 105.3021240234375 | -3.25 | 0.68179327336871909 | 0.076297157689144732 | 0.75 | -3.5398201437373751 | 0.0029706284773525033 | 1.8229790322483097e-24 | 4.6184756093705365e-13 | PASS |
| 48 | vascular_length | 516.375 | 557.51202769480824 | -11.6875 | 1.5803680032336627 | 0.079665025794835698 | 0.75 | 10.620622989835244 | 2.2500667964803482e-08 | 7.0984742069798786e-30 | 1.2866063866318144e-18 | PASS |
| 48 | perfused_volume | 250.21955695785354 | 276.85969750324449 | -4.422932266954918 | 2.5690258202060288 | 0.10646705984647788 | 0.75 | 17.126812052972241 | 2.94634644604435e-11 | 2.2420379513499188e-31 | 1.8336890864437284e-19 | PASS |
| 48 | lesion_perfused_fraction | 0.41740702086600523 | 0.46172580229214238 | -0.0069021640992518833 | 2.6110925335295416 | 0.044318781426137113 | 0.14999999999999999 | 17.245313849840532 | 2.6693724456472291e-11 | 4.3582536421269156e-21 | 3.8728107352972257e-17 | PASS |
| 48 | vascular_branches | 481.1875 | 550.89060094484125 | 6.1875 | 0.96628876074687375 | 0.144856424875628 | 0.75 | 5.4614466805701749 | 6.5634621630142003e-05 | 1.5499391492534712e-23 | 5.3390475520226672e-10 | PASS |
| 48 | vascular_anastomoses | 587.8125 | 590.8897929277806 | -0.5625 | 0.041737517878473955 | 0.005235160749015394 | 0.75 | 0.25072373411368992 | 0.80542959207420184 | 2.325288900177335e-24 | 1.0191408352808234e-12 | PASS |
| 48 | vascular_roots | 138.625 | 114.5220947265625 | -5.4375 | 1.8024800075268563 | 0.17387127338818756 | 0.75 | -9.0859152164387122 | 1.7367146994897871e-07 | 5.041809053599588e-24 | 1.5752405395895139e-14 | PASS |
| 72 | vascular_length | 528.125 | 596.95613785108571 | -6.75 | 3.2786026718935997 | 0.13033114859377176 | 0.75 | 17.419000125333561 | 2.3123332974958627e-11 | 4.8532036522028329e-30 | 5.0862880259775603e-18 | PASS |
| 72 | perfused_volume | 255.38261680919607 | 299.72456785486656 | -1.5553122014586016 | 5.2842310983542466 | 0.17362948034477868 | 0.75 | 26.009855890744426 | 6.80835967933087e-14 | 4.5245410928199555e-31 | 3.1184102173508266e-18 | PASS |
| 72 | lesion_perfused_fraction | 0.42602633944702595 | 0.49834746032752852 | -0.0021828725078417745 | 5.0408219515858121 | 0.072321120880502562 | 0.14999999999999999 | 25.711402650108852 | 8.0662664610450413e-14 | 2.2456587882757677e-21 | 1.4097728664114604e-14 | PASS |
| 72 | vascular_branches | 482.5 | 553.15036696282493 | 7.3125 | 0.98584773268076265 | 0.14642563101103606 | 0.75 | 5.5036988215877907 | 6.063208407378082e-05 | 1.6155909247625661e-23 | 5.7779096628155846e-10 | PASS |
| 72 | vascular_anastomoses | 618.375 | 641.2591499173335 | 0.375 | 0.32566848405641935 | 0.037006913147092764 | 0.75 | 1.8440588309095309 | 0.085016973252516351 | 1.4805700506140486e-24 | 1.0323956246007184e-12 | PASS |
| 72 | vascular_roots | 139.25 | 114.5220947265625 | -4.6875 | 1.8524594539555055 | 0.17757921201750448 | 0.75 | -9.3748388124815314 | 1.1603065546943496e-07 | 4.1641167057397811e-24 | 1.3042752656061224e-14 | PASS |
| 96 | vascular_length | 529.5 | 621.69958077175113 | -6.125 | 4.5177421865965117 | 0.17412574272285397 | 0.75 | 22.986214416512194 | 4.1639774454514769e-13 | 4.6800151095096541e-30 | 1.9696736546416587e-17 | PASS |
| 96 | perfused_volume | 256.08558524338457 | 316.79429369929881 | -1.3424184393529508 | 7.33819990166569 | 0.23706413774994983 | 0.75 | 34.638462756351657 | 9.9055180025722325e-16 | 6.4155930062374187e-31 | 2.631321444876912e-17 | PASS |
| 96 | lesion_perfused_fraction | 0.42715638396456895 | 0.52242772708287044 | -0.0017759324199288704 | 6.6573496596981911 | 0.095271343118301488 | 0.14999999999999999 | 32.838749409374486 | 2.1837209261919244e-15 | 8.2001525908572104e-22 | 3.6740260455349311e-12 | PASS |
| 96 | vascular_branches | 482.5 | 553.20519467863221 | 7.3125 | 0.98661279281619219 | 0.14653926358265748 | 0.75 | 5.5079103537503586 | 6.0155622063913377e-05 | 1.6148857781647738e-23 | 5.793369840628025e-10 | PASS |
| 96 | vascular_anastomoses | 621.4375 | 656.31083442687418 | 2.625 | 0.49857856749531737 | 0.056117203140901809 | 0.75 | 2.7980154565014508 | 0.013513724029793218 | 1.2741635651512747e-24 | 1.4972629853195224e-12 | PASS |
| 96 | vascular_roots | 139.25 | 114.5220947265625 | -4.625 | 1.85735195156991 | 0.17757921201750448 | 0.75 | -9.3748388124815314 | 1.1603065546943496e-07 | 4.1641167057397811e-24 | 1.3042752656061224e-14 | PASS |
| 120 | vascular_length | 529.875 | 639.198522691937 | -6.25 | 5.3762054948665972 | 0.20631945778143343 | 0.75 | 27.520237442530622 | 2.9678164248017997e-14 | 4.1922151532239061e-30 | 4.1540866424378927e-17 | PASS |
| 120 | perfused_volume | 256.24645944835868 | 330.08977061624165 | -1.3822213831228627 | 9.1354495846251531 | 0.28817300081667913 | 0.75 | 42.729656480991082 | 4.3784640465609621e-17 | 7.1008265681172565e-31 | 1.0044281581340728e-16 | PASS |
| 120 | lesion_perfused_fraction | 0.42742288846285914 | 0.53789708140407833 | -0.0018385290064696046 | 7.8734470148090647 | 0.11047419294121917 | 0.14999999999999999 | 38.459603156628035 | 2.09575832659777e-16 | 2.8714337441029472e-22 | 3.2596296036677254e-10 | PASS |
| 120 | vascular_branches | 482.5 | 553.206907468802 | 7.3125 | 0.98663669290868961 | 0.14654281340684347 | 0.75 | 5.5080420182571714 | 6.0140789293394366e-05 | 1.614860494639758e-23 | 5.7938528706737084e-10 | PASS |
| 120 | vascular_anastomoses | 621.625 | 661.85046999595659 | 2.75 | 0.57373554252538084 | 0.064710187003348602 | 0.75 | 3.2136128251937186 | 0.0057999078526985897 | 1.2098483423538933e-24 | 1.8997598792181036e-12 | PASS |
| 120 | vascular_roots | 139.25 | 114.5220947265625 | -4.625 | 1.85735195156991 | 0.17757921201750448 | 0.75 | -9.3748388124815314 | 1.1603065546943496e-07 | 4.1641167057397811e-24 | 1.3042752656061224e-14 | PASS |
| 144 | vascular_length | 529.9375 | 652.36817124660331 | -6.25 | 5.9909678124062635 | 0.23102851043114195 | 0.75 | 30.802374901924562 | 5.6298063613898567e-15 | 4.3087954230240915e-30 | 8.3093102572956032e-17 | PASS |
| 144 | perfused_volume | 256.26099345822962 | 340.72943846595183 | -1.378661732865913 | 10.425987745898453 | 0.32961881505189178 | 0.75 | 48.773867999146603 | 6.0991347344233873e-18 | 8.4885778198370365e-31 | 3.9712278014959813e-16 | PASS |
| 144 | lesion_perfused_fraction | 0.42744691161967058 | 0.54778375030475535 | -0.0018323457189941134 | 8.570868919109488 | 0.12033683868508473 | 0.14999999999999999 | 41.701138370495876 | 6.292118856614898e-17 | 1.7630419979063066e-22 | 1.7366137214171909e-08 | PASS |
| 144 | vascular_branches | 482.5 | 553.20697399379617 | 7.3125 | 0.9866376211914073 | 0.14654295128247902 | 0.75 | 5.5080471355363949 | 6.0140212878836991e-05 | 1.6148593699272721e-23 | 5.7938716138137153e-10 | PASS |
| 144 | vascular_anastomoses | 621.6875 | 664.16126349301635 | 2.75 | 0.6055104502148152 | 0.068320118215367598 | 0.75 | 3.3976512941270167 | 0.0039769030386126525 | 1.1484632534007238e-24 | 2.0064449011951331e-12 | PASS |
| 144 | vascular_roots | 139.25 | 114.5220947265625 | -4.625 | 1.85735195156991 | 0.17757921201750448 | 0.75 | -9.3748388124815314 | 1.1603065546943496e-07 | 4.1641167057397811e-24 | 1.3042752656061224e-14 | PASS |
| 168 | vascular_length | 530 | 662.75472012532998 | -6.3125 | 6.4533252743275389 | 0.25048060401005667 | 0.75 | 33.345289118311356 | 1.7407939997170611e-15 | 4.319893122051547e-30 | 1.4936304681927387e-16 | PASS |
| 168 | perfused_volume | 256.30050099315639 | 349.47241552304587 | -1.4181692677926971 | 11.385462503552883 | 0.36352607259389358 | 0.75 | 53.540789379546581 | 1.5173114220883127e-18 | 9.2790335484236622e-31 | 1.4301915170938175e-15 | PASS |
| 168 | lesion_perfused_fraction | 0.42751221333029338 | 0.55420589157645661 | -0.0018976474296168964 | 8.9677115398752765 | 0.12669367824616323 | 0.14999999999999999 | 43.610010297207452 | 3.2319791571063e-17 | 1.3765551979735451e-22 | 4.1595766225275137e-07 | PASS |
| 168 | vascular_branches | 482.5 | 553.20697714116159 | 7.3125 | 0.98663766510941142 | 0.14654295780551624 | 0.75 | 5.5080473777819572 | 6.0140185592239084e-05 | 1.6148593096825189e-23 | 5.7938724999136801e-10 | PASS |
| 168 | vascular_anastomoses | 621.75 | 665.24868803422521 | 2.6875 | 0.61944164603353957 | 0.06996170170361922 | 0.75 | 3.4844630025969945 | 0.003327960397220644 | 1.1088345958752569e-24 | 2.0326849667029615e-12 | PASS |
| 168 | vascular_roots | 139.25 | 114.5220947265625 | -4.625 | 1.85735195156991 | 0.17757921201750448 | 0.75 | -9.3748388124815314 | 1.1603065546943496e-07 | 4.1641167057397811e-24 | 1.3042752656061224e-14 | PASS |
| 192 | vascular_length | 530 | 671.2504818118415 | -6.3125 | 6.8663118224854376 | 0.26651034304121035 | 0.75 | 35.498585211268413 | 6.8849919341279519e-16 | 4.3070124664664933e-30 | 2.3763062437407226e-16 | PASS |
| 192 | perfused_volume | 256.30050099315639 | 356.82919907462832 | -1.4181692677926971 | 12.284449968776933 | 0.39222981497081127 | 0.75 | 57.684948139642366 | 4.9844819942399615e-19 | 9.8282227073802734e-31 | 4.383823187916353e-15 | PASS |
| 192 | lesion_perfused_fraction | 0.42751221333029338 | 0.55848758831187562 | -0.0018976474296168964 | 9.2707812885478109 | 0.13097537498158227 | 0.14999999999999999 | 44.982564174832412 | 2.0372312555440345e-17 | 1.1312003597572477e-22 | 4.7289265773061439e-06 | PASS |
| 192 | vascular_branches | 482.5 | 553.20697731394284 | 7.3125 | 0.98663766752038251 | 0.14654295816361207 | 0.75 | 5.5080473910847179 | 6.0140184093813202e-05 | 1.6148593061285572e-23 | 5.7938725485427732e-10 | PASS |
| 192 | vascular_anastomoses | 621.75 | 665.83117435974066 | 2.6875 | 0.6277365235246809 | 0.07089855144308907 | 0.75 | 3.5309631993362234 | 0.0030251028944132115 | 1.0989558911441757e-24 | 2.0740275542145248e-12 | PASS |
| 192 | vascular_roots | 139.25 | 114.5220947265625 | -4.625 | 1.85735195156991 | 0.17757921201750448 | 0.75 | -9.3748388124815314 | 1.1603065546943496e-07 | 4.1641167057397811e-24 | 1.3042752656061224e-14 | PASS |
| 216 | vascular_length | 530 | 678.39855208033634 | -6.3125 | 7.2137858895679345 | 0.27999726807610631 | 0.75 | 37.306732659889491 | 3.2936232623043424e-16 | 4.245030250205059e-30 | 3.5627128723528957e-16 | PASS |
| 216 | perfused_volume | 256.30050099315639 | 363.14519523266404 | -1.4181692677926971 | 13.056254839297553 | 0.41687274829931209 | 0.75 | 61.223382882489076 | 2.0480590329895016e-19 | 1.0081400634927778e-30 | 1.2370489259513649e-14 | PASS |
| 216 | lesion_perfused_fraction | 0.42751221333029338 | 0.56142663980392937 | -0.0018976474296168964 | 9.4788150779715199 | 0.13391442647363594 | 0.14999999999999999 | 45.912109738458099 | 1.5021008269171041e-17 | 9.9347365626208011e-23 | 2.9687563299179185e-05 | PASS |
| 216 | vascular_branches | 482.5 | 553.20697732051917 | 7.3125 | 0.98663766761214811 | 0.14654295817724181 | 0.75 | 5.5080473915885673 | 6.0140184037059548e-05 | 1.6148593060790492e-23 | 5.7938725504091036e-10 | PASS |
| 216 | vascular_anastomoses | 621.75 | 666.18476960162627 | 2.6875 | 0.63277188501629489 | 0.071467261120428255 | 0.75 | 3.5591914879782527 | 0.0028548894872594522 | 1.092662966728046e-24 | 2.0995652399490627e-12 | PASS |
| 216 | vascular_roots | 139.25 | 114.5220947265625 | -4.625 | 1.85735195156991 | 0.17757921201750448 | 0.75 | -9.3748388124815314 | 1.1603065546943496e-07 | 4.1641167057397811e-24 | 1.3042752656061224e-14 | PASS |
| 240 | vascular_length | 530 | 684.54572450997489 | -6.3125 | 7.5126054205675432 | 0.29159570662259404 | 0.75 | 38.859399941119975 | 1.7972229141980815e-16 | 4.1542381837643746e-30 | 5.1018897898563699e-16 | PASS |
| 240 | perfused_volume | 256.30050099315639 | 368.65860492795275 | -1.4181692677926971 | 13.729984896999621 | 0.43838425402764436 | 0.75 | 64.299041901454729 | 9.8441004187909862e-20 | 1.0128501730834142e-30 | 3.2652664822662161e-14 | PASS |
| 240 | lesion_perfused_fraction | 0.42751221333029338 | 0.56350292394615076 | -0.0018976474296168964 | 9.6257799267312336 | 0.13599071061585738 | 0.14999999999999999 | 46.561820617375872 | 1.218313049039064e-17 | 9.0872886388442132e-23 | 0.00011771753757838816 | PASS |
| 240 | vascular_branches | 482.5 | 553.20697732051917 | 7.3125 | 0.98663766761214811 | 0.14654295817724181 | 0.75 | 5.5080473915885673 | 6.0140184037059548e-05 | 1.6148593060790492e-23 | 5.7938725504091036e-10 | PASS |
| 240 | vascular_anastomoses | 621.75 | 666.4229066070875 | 2.6875 | 0.63616306726363814 | 0.071850271985665418 | 0.75 | 3.5782038799228149 | 0.0027456950586177264 | 1.0882131800779266e-24 | 2.1169538591427752e-12 | PASS |
| 240 | vascular_roots | 139.25 | 114.5220947265625 | -4.625 | 1.85735195156991 | 0.17757921201750448 | 0.75 | -9.3748388124815314 | 1.1603065546943496e-07 | 4.1641167057397811e-24 | 1.3042752656061224e-14 | PASS |

## 2000-square field benchmark

The advance-only benchmark runs 24 hours with 0.25-hour steps and a 128-square
frozen consumer patch. Allocation is excluded from advance time; process peak
RSS includes inputs and reused work buffers. Three separate processes per
thread count use identical initial conditions. This is not a whole-model
production benchmark. The raw records follow.

```json
{
  "schema_version": 1,
  "system": "Darwin",
  "release": "27.0.0",
  "machine": "arm64",
  "records": [
    {
      "grid_shape": [
        2000,
        2000,
        1
      ],
      "threads": 1,
      "hours": 24,
      "step_hours": 0.25,
      "updates": 96,
      "source_edge": 128,
      "wall_seconds": 1.41893525,
      "peak_resident_bytes": 643219456,
      "allocated_field_bytes": 544008192,
      "active_voxels": 51076,
      "vascular_length": 1611.0819832119057,
      "vascular_branches": 627.651060860784,
      "vascular_anastomoses": 74.18064820707093,
      "state_checksum": 10645225464068239034,
      "repetition": 1
    },
    {
      "grid_shape": [
        2000,
        2000,
        1
      ],
      "threads": 8,
      "hours": 24,
      "step_hours": 0.25,
      "updates": 96,
      "source_edge": 128,
      "wall_seconds": 0.5675905,
      "peak_resident_bytes": 643612672,
      "allocated_field_bytes": 544008192,
      "active_voxels": 51076,
      "vascular_length": 1611.0819832119057,
      "vascular_branches": 627.651060860784,
      "vascular_anastomoses": 74.18064820707093,
      "state_checksum": 10645225464068239034,
      "repetition": 1
    },
    {
      "grid_shape": [
        2000,
        2000,
        1
      ],
      "threads": 1,
      "hours": 24,
      "step_hours": 0.25,
      "updates": 96,
      "source_edge": 128,
      "wall_seconds": 1.302991125,
      "peak_resident_bytes": 643203072,
      "allocated_field_bytes": 544008192,
      "active_voxels": 51076,
      "vascular_length": 1611.0819832119057,
      "vascular_branches": 627.651060860784,
      "vascular_anastomoses": 74.18064820707093,
      "state_checksum": 10645225464068239034,
      "repetition": 2
    },
    {
      "grid_shape": [
        2000,
        2000,
        1
      ],
      "threads": 8,
      "hours": 24,
      "step_hours": 0.25,
      "updates": 96,
      "source_edge": 128,
      "wall_seconds": 0.474532834,
      "peak_resident_bytes": 643563520,
      "allocated_field_bytes": 544008192,
      "active_voxels": 51076,
      "vascular_length": 1611.0819832119057,
      "vascular_branches": 627.651060860784,
      "vascular_anastomoses": 74.18064820707093,
      "state_checksum": 10645225464068239034,
      "repetition": 2
    },
    {
      "grid_shape": [
        2000,
        2000,
        1
      ],
      "threads": 1,
      "hours": 24,
      "step_hours": 0.25,
      "updates": 96,
      "source_edge": 128,
      "wall_seconds": 1.215423417,
      "peak_resident_bytes": 643219456,
      "allocated_field_bytes": 544008192,
      "active_voxels": 51076,
      "vascular_length": 1611.0819832119057,
      "vascular_branches": 627.651060860784,
      "vascular_anastomoses": 74.18064820707093,
      "state_checksum": 10645225464068239034,
      "repetition": 3
    },
    {
      "grid_shape": [
        2000,
        2000,
        1
      ],
      "threads": 8,
      "hours": 24,
      "step_hours": 0.25,
      "updates": 96,
      "source_edge": 128,
      "wall_seconds": 0.491949792,
      "peak_resident_bytes": 643563520,
      "allocated_field_bytes": 544008192,
      "active_voxels": 51076,
      "vascular_length": 1611.0819832119057,
      "vascular_branches": 627.651060860784,
      "vascular_anastomoses": 74.18064820707093,
      "state_checksum": 10645225464068239034,
      "repetition": 3
    }
  ],
  "deterministic": true,
  "median_wall_seconds": {
    "1": 1.302991125,
    "8": 0.491949792
  },
  "speedup_1_to_8": 2.648626234199119,
  "timing_scope": "field advance only; peak resident memory includes all buffers and inputs"
}
```

## Remaining limitations

- The ensemble is a frozen-lesion mechanism experiment. A growing invasive
  lesion and an angiogenic hybrid ensemble remain unvalidated.
- Roots are perfused host-supply inlets. Connectivity to a remote vascular
  network is not solved. Connection counts are local tip termination events,
  including crowding, not reconstructed unique geometric junction counts.
- PDE independent-arrival edge closure loses repeated-path correlations.
  At 240 hours mean centerline is 530 versus 684.5457245099749 and mean
  perfused volume is 256.3005009931564 versus 368.65860492795275.
- Fractional field deletion and indivisible ABM footprint deletion can differ.
  These losses are diagnosed; there is no conservative relocation claim.
- The existing long phase B cases remain opt-in and unexecuted. Three-dimensional
  production performance, front-based hybrid classification and operator strategy
  refactoring are separate later work.
- Local Apple clang results do not certify remote Linux/macOS CI completion.
