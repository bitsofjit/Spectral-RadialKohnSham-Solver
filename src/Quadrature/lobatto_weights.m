function w = lobatto_weights(n)

w = ones(1,n);
switch n
   case(2)
      w(1)=1.0000000000000000000000000000000;
      w(2)=1.0000000000000000000000000000000;
   case(3)
      w(1)=0.33333333333333333333333333333333;
      w(2)=1.3333333333333333333333333333333;
      w(3)=0.33333333333333333333333333333333;
   case(4)
      w(1)=0.16666666666666666666666666666667;
      w(2)=0.83333333333333333333333333333333;
      w(3)=0.83333333333333333333333333333333;
      w(4)=0.16666666666666666666666666666667;
   case(5)
      w(1)=0.10000000000000000000000000000000;
      w(2)=0.54444444444444444444444444444444;
      w(3)=0.71111111111111111111111111111111;
      w(4)=0.54444444444444444444444444444444;
      w(5)=0.10000000000000000000000000000000;
   case(6)
      w(1)=0.066666666666666666666666666666667;
      w(2)=0.37847495629784698031661280821202;
      w(3)=0.55485837703548635301672052512131;
      w(4)=0.55485837703548635301672052512131;
      w(5)=0.37847495629784698031661280821202;
      w(6)=0.066666666666666666666666666666667;
   case(7)
      w(1)=0.047619047619047619047619047619048;
      w(2)=0.27682604736156594801070040629007;
      w(3)=0.43174538120986262341787102228136;
      w(4)=0.48761904761904761904761904761905;
      w(5)=0.43174538120986262341787102228136;
      w(6)=0.27682604736156594801070040629007;
      w(7)=0.047619047619047619047619047619048;
   case(8)
      w(1)=0.035714285714285714285714285714286;
      w(2)=0.21070422714350603938299206577576;
      w(3)=0.34112269248350436476424067710775;
      w(4)=0.41245879465870388156705297140221;
      w(5)=0.41245879465870388156705297140221;
      w(6)=0.34112269248350436476424067710775;
      w(7)=0.21070422714350603938299206577576;
      w(8)=0.035714285714285714285714285714286;
   case(9)
      w(1)=0.027777777777777777777777777777778;
      w(2)=0.16549536156080552504633972002921;
      w(3)=0.27453871250016173528070561857937;
      w(4)=0.34642851097304634511513153213972;
      w(5)=0.37151927437641723356009070294785;
      w(6)=0.34642851097304634511513153213972;
      w(7)=0.27453871250016173528070561857937;
      w(8)=0.16549536156080552504633972002921;
      w(9)=0.027777777777777777777777777777778;
   case(10)
      w(1)=0.022222222222222222222222222222222;
      w(2)=0.13330599085107011112622717075539;
      w(3)=0.22488934206312645211945782173105;
      w(4)=0.29204268367968375787558225737444;
      w(5)=0.32753976118389745665651052791689;
      w(6)=0.32753976118389745665651052791689;
      w(7)=0.29204268367968375787558225737444;
      w(8)=0.22488934206312645211945782173105;
      w(9)=0.13330599085107011112622717075539;
      w(10)=0.022222222222222222222222222222222;
   case(11)
      w(1)=0.018181818181818181818181818181818;
      w(2)=0.10961227326699486446140344958035;
      w(3)=0.18716988178030520410814152189943;
      w(4)=0.24804810426402831404008486642187;
      w(5)=0.28687912477900808867922240333154;
      w(6)=0.30021759545569069378593188116998;
      w(7)=0.28687912477900808867922240333154;
      w(8)=0.24804810426402831404008486642187;
      w(9)=0.18716988178030520410814152189943;
      w(10)=0.10961227326699486446140344958035;
      w(11)=0.018181818181818181818181818181818;
   case(12)
      w(1)=0.015151515151515151515151515151515;
      w(2)=0.091684517413196130668342594134079;
      w(3)=0.15797470556437011516467106270034;
      w(4)=0.21250841776102114535830207736687;
      w(5)=0.25127560319920128029324441214760;
      w(6)=0.27140524091069617700028833849960;
      w(7)=0.27140524091069617700028833849960;
      w(8)=0.25127560319920128029324441214760;
      w(9)=0.21250841776102114535830207736687;
      w(10)=0.15797470556437011516467106270034;
      w(11)=0.091684517413196130668342594134079;
      w(12)=0.015151515151515151515151515151515;
   case(13)
      w(1)=0.012820512820512820512820512820513;
      w(2)=0.077801686746818927793588988333134;
      w(3)=0.13498192668960834911991476258937;
      w(4)=0.18364686520355009200749425874681;
      w(5)=0.22076779356611008608553400837940;
      w(6)=0.24401579030667635645857814836016;
      w(7)=0.25193084933344673604413864154124;
      w(8)=0.24401579030667635645857814836016;
      w(9)=0.22076779356611008608553400837940;
      w(10)=0.18364686520355009200749425874681;
      w(11)=0.13498192668960834911991476258937;
      w(12)=0.077801686746818927793588988333134;
      w(13)=0.012820512820512820512820512820513;
   case(14)
      w(1)=0.010989010989010989010989010989011;
      w(2)=0.066837284497681284634070660746053;
      w(3)=0.11658665589871165154099667065465;
      w(4)=0.16002185176295214241282099798759;
      w(5)=0.19482614937341611864033177837588;
      w(6)=0.21912625300977075487116252395417;
      w(7)=0.23161279446845705888962835729264;
      w(8)=0.23161279446845705888962835729264;
      w(9)=0.21912625300977075487116252395417;
      w(10)=0.19482614937341611864033177837588;
      w(11)=0.16002185176295214241282099798759;
      w(12)=0.11658665589871165154099667065465;
      w(13)=0.066837284497681284634070660746053;
      w(14)=0.010989010989010989010989010989011;
   case(15)
      w(1)=0.0095238095238095238095238095238095;
      w(2)=0.058029893028601249096880584025282;
      w(3)=0.10166007032571806760366617078880;
      w(4)=0.14051169980242810946044680564367;
      w(5)=0.17278964725360094905207709940835;
      w(6)=0.19698723596461335609250034650741;
      w(7)=0.21197358592682092012743007697722;
      w(8)=0.21704811634881564951495021425091;
      w(9)=0.21197358592682092012743007697722;
      w(10)=0.19698723596461335609250034650741;
      w(11)=0.17278964725360094905207709940835;
      w(12)=0.14051169980242810946044680564367;
      w(13)=0.10166007032571806760366617078880;
      w(14)=0.058029893028601249096880584025282;
      w(15)=0.0095238095238095238095238095238095;
   case(16)
      w(1)=0.0083333333333333333333333333333333;
      w(2)=0.050850361005919905403244919565455;
      w(3)=0.089393697325930800991052080166084;
      w(4)=0.12425538213251409834953633265731;
      w(5)=0.15402698080716428081564494048499;
      w(6)=0.17749191339170412530107566952836;
      w(7)=0.19369002382520358431691359885352;
      w(8)=0.20195830817822987148919912541094;
      w(9)=0.20195830817822987148919912541094;
      w(10)=0.19369002382520358431691359885352;
      w(11)=0.17749191339170412530107566952836;
      w(12)=0.15402698080716428081564494048499;
      w(13)=0.12425538213251409834953633265731;
      w(14)=0.089393697325930800991052080166084;
      w(15)=0.050850361005919905403244919565455;
      w(16)=0.0083333333333333333333333333333333;
   case(17)
      w(1)=0.0073529411764705882352941176470588;
      w(2)=0.044921940543254209647400954623212;
      w(3)=0.079198270503687119190264429952835;
      w(4)=0.11059290900702816137577270522008;
      w(5)=0.13798774620192655905620157495403;
      w(6)=0.16039466199762153951632836586475;
      w(7)=0.17700425351565787043694574536329;
      w(8)=0.18721633967761923589208848286062;
      w(9)=0.19066187475346943329940724702825;
      w(10)=0.18721633967761923589208848286062;
      w(11)=0.17700425351565787043694574536329;
      w(12)=0.16039466199762153951632836586475;
      w(13)=0.13798774620192655905620157495403;
      w(14)=0.11059290900702816137577270522008;
      w(15)=0.079198270503687119190264429952835;
      w(16)=0.044921940543254209647400954623212;
      w(17)=0.0073529411764705882352941176470588;
   case(18)
      w(1)=0.0065359477124183006535947712418301;
      w(2)=0.039970628810914066137599176410101;
      w(3)=0.070637166885633664999222960167786;
      w(4)=0.099016271717502802394423605318672;
      w(5)=0.12421053313296710026339635889675;
      w(6)=0.14541196157380226798300321049443;
      w(7)=0.16193951723760248926432670670023;
      w(8)=0.17326210948945622601061440382668;
      w(9)=0.17901586343970308229381880694353;
      w(10)=0.17901586343970308229381880694353;
      w(11)=0.17326210948945622601061440382668;
      w(12)=0.16193951723760248926432670670023;
      w(13)=0.14541196157380226798300321049443;
      w(14)=0.12421053313296710026339635889675;
      w(15)=0.099016271717502802394423605318672;
      w(16)=0.070637166885633664999222960167786;
      w(17)=0.039970628810914066137599176410101;
      w(18)=0.0065359477124183006535947712418301;
   case(19)
      w(1)=0.0058479532163742690058479532163743;
      w(2)=0.035793365186176477115425569035122;
      w(3)=0.063381891762629736851695690418317;
      w(4)=0.089131757099207084448008790556153;
      w(5)=0.11231534147730504407091001546378;
      w(6)=0.13226728044875077692604673390973;
      w(7)=0.14841394259593888500968064366841;
      w(8)=0.16029092404406124197991096818359;
      w(9)=0.16755658452714286727013727774026;
      w(10)=0.17000191928482723464467271561652;
      w(11)=0.16755658452714286727013727774026;
      w(12)=0.16029092404406124197991096818359;
      w(13)=0.14841394259593888500968064366841;
      w(14)=0.13226728044875077692604673390973;
      w(15)=0.11231534147730504407091001546378;
      w(16)=0.089131757099207084448008790556153;
      w(17)=0.063381891762629736851695690418317;
      w(18)=0.035793365186176477115425569035122;
      w(19)=0.0058479532163742690058479532163743;
   case(20)
      w(1)=0.0052631578947368421052631578947368;
      w(2)=0.032237123188488941491605028117294;
      w(3)=0.057181802127566826004753627173243;
      w(4)=0.080631763996119603144776846113721;
      w(5)=0.10199149969945081568378120573289;
      w(6)=0.12070922762867472509942970500239;
      w(7)=0.13630048235872418448978079298903;
      w(8)=0.14836155407091682581471301373397;
      w(9)=0.15658010264747548715816989679364;
      w(10)=0.16074328638784574900772672644908;
      w(11)=0.16074328638784574900772672644908;
      w(12)=0.15658010264747548715816989679364;
      w(13)=0.14836155407091682581471301373397;
      w(14)=0.13630048235872418448978079298903;
      w(15)=0.12070922762867472509942970500239;
      w(16)=0.10199149969945081568378120573289;
      w(17)=0.080631763996119603144776846113721;
      w(18)=0.057181802127566826004753627173243;
      w(19)=0.032237123188488941491605028117294;
      w(20)=0.0052631578947368421052631578947368;
   case(21)
      w(1)=0.0047619047619047619047619047619048;
      w(2)=0.029184840098505458609458543613171;
      w(3)=0.051843169000849625072722971852830;
      w(4)=0.073273918185074144252547861041894;
      w(5)=0.092985467957886065301137664149214;
      w(6)=0.11051708321912333526700048678439;
      w(7)=0.12545812119086894801515753570800;
      w(8)=0.13745846286004134358089961741515;
      w(9)=0.14623686244797745926727053063439;
      w(10)=0.15158757511168138445325068150529;
      w(11)=0.15338519033217494855158440506754;
      w(12)=0.15158757511168138445325068150529;
      w(13)=0.14623686244797745926727053063439;
      w(14)=0.13745846286004134358089961741515;
      w(15)=0.12545812119086894801515753570800;
      w(16)=0.11051708321912333526700048678439;
      w(17)=0.092985467957886065301137664149214;
      w(18)=0.073273918185074144252547861041894;
      w(19)=0.051843169000849625072722971852830;
      w(20)=0.029184840098505458609458543613171;
      w(21)=0.0047619047619047619047619047619048;
   case(22)
      w(1)=0.0043290043290043290043290043290043;
      w(2)=0.026545747682501757911627904520543;
      w(3)=0.047214465293740752123775734864792;
      w(4)=0.066865605864553076012404194157097;
      w(5)=0.085090060391838447815711236095748;
      w(6)=0.10150057480164767437243730374960;
      w(7)=0.11574764465393906659003636772146;
      w(8)=0.12752769665343027553084445930883;
      w(9)=0.13658968861374142668617736220617;
      w(10)=0.14274049227136140033623599356679;
      w(11)=0.14584901944424179361642043947997;
      w(12)=0.14584901944424179361642043947997;
      w(13)=0.14274049227136140033623599356679;
      w(14)=0.13658968861374142668617736220617;
      w(15)=0.12752769665343027553084445930883;
      w(16)=0.11574764465393906659003636772146;
      w(17)=0.10150057480164767437243730374960;
      w(18)=0.085090060391838447815711236095748;
      w(19)=0.066865605864553076012404194157097;
      w(20)=0.047214465293740752123775734864792;
      w(21)=0.026545747682501757911627904520543;
      w(22)=0.0043290043290043290043290043290043;
   case(23)
      w(1)=0.0039525691699604743083003952569170;
      w(2)=0.024248600771531736517399658937097;
      w(3)=0.043175871170241834748876465612042;
      w(4)=0.061252477129554206381382847440355;
      w(5)=0.078135449475569989741934255347965;
      w(6)=0.093497246163512341833500706906697;
      w(7)=0.10703910172433651153518362791547;
      w(8)=0.11849751066274913130212600472426;
      w(9)=0.12764947470175887663614855305567;
      w(10)=0.13431687263860381990156489770071;
      w(11)=0.13836993638580739452350273386294;
      w(12)=0.13972978001274736514015970647975;
      w(13)=0.13836993638580739452350273386294;
      w(14)=0.13431687263860381990156489770071;
      w(15)=0.12764947470175887663614855305567;
      w(16)=0.11849751066274913130212600472426;
      w(17)=0.10703910172433651153518362791547;
      w(18)=0.093497246163512341833500706906697;
      w(19)=0.078135449475569989741934255347965;
      w(20)=0.061252477129554206381382847440355;
      w(21)=0.043175871170241834748876465612042;
      w(22)=0.024248600771531736517399658937097;
      w(23)=0.0039525691699604743083003952569170;
   case(24)
      w(1)=0.0036231884057971014492753623188406;
      w(2)=0.022236853464711208992960434817299;
      w(3)=0.039631681333467809469262163363472;
      w(4)=0.056309848724646199020948339678118;
      w(5)=0.071981862055293982215639060668680;
      w(6)=0.086369029967929068216562819306583;
      w(7)=0.099214827684083587414149770995712;
      w(8)=0.11029008689296860411042412302760;
      w(9)=0.11939719370249131903229388835072;
      w(10)=0.12637364202802080012724809069137;
      w(11)=0.13109494187360394235445022581653;
      w(12)=0.13347684386698637759678572096508;
      w(13)=0.13347684386698637759678572096508;
      w(14)=0.13109494187360394235445022581653;
      w(15)=0.12637364202802080012724809069137;
      w(16)=0.11939719370249131903229388835072;
      w(17)=0.11029008689296860411042412302760;
      w(18)=0.099214827684083587414149770995712;
      w(19)=0.086369029967929068216562819306583;
      w(20)=0.071981862055293982215639060668680;
      w(21)=0.056309848724646199020948339678118;
      w(22)=0.039631681333467809469262163363472;
      w(23)=0.022236853464711208992960434817299;
      w(24)=0.0036231884057971014492753623188406;
   case(25)
      w(1)=0.0033333333333333333333333333333333;
      w(2)=0.020465168932974385308542471819163;
      w(3)=0.036504738794271372032382988755110;
      w(4)=0.051936228368491474643333889713489;
      w(5)=0.066513728675312784693869993313599;
      w(6)=0.079998774836292981801626436138836;
      w(7)=0.092170139910620421912689622714049;
      w(8)=0.10282803034795783082750364490713;
      w(9)=0.11179746626832088815624423243104;
      w(10)=0.11893117940681182540944424463419;
      w(11)=0.12411203893795029069521531243930;
      w(12)=0.12725497753833144701711441673876;
      w(13)=0.12830838929866192833739882612401;
      w(14)=0.12725497753833144701711441673876;
      w(15)=0.12411203893795029069521531243930;
      w(16)=0.11893117940681182540944424463419;
      w(17)=0.11179746626832088815624423243104;
      w(18)=0.10282803034795783082750364490713;
      w(19)=0.092170139910620421912689622714049;
      w(20)=0.079998774836292981801626436138836;
      w(21)=0.066513728675312784693869993313599;
      w(22)=0.051936228368491474643333889713489;
      w(23)=0.036504738794271372032382988755110;
      w(24)=0.020465168932974385308542471819163;
      w(25)=0.0033333333333333333333333333333333;
   case(26)
      w(1)=0.0030769230769230769230769230769231;
      w(2)=0.018896858024263465581352047676267;
      w(3)=0.033732303685955999377522748582756;
      w(4)=0.048048399081180627315979312780662;
      w(5)=0.061635025142547402782180054104279;
      w(6)=0.074287050122291137316084082986249;
      w(7)=0.085812863980004362187481252359068;
      w(8)=0.096037802353901310803697396370862;
      w(9)=0.10480688623073705298990583689136;
      w(10)=0.11198719411986033530021087987989;
      w(11)=0.11746988409380900707669562029648;
      w(12)=0.12117184628844334604428669570423;
      w(13)=0.12303696380008287630152714929097;
      w(14)=0.12303696380008287630152714929097;
      w(15)=0.12117184628844334604428669570423;
      w(16)=0.11746988409380900707669562029648;
      w(17)=0.11198719411986033530021087987989;
      w(18)=0.10480688623073705298990583689136;
      w(19)=0.096037802353901310803697396370862;
      w(20)=0.085812863980004362187481252359068;
      w(21)=0.074287050122291137316084082986249;
      w(22)=0.061635025142547402782180054104279;
      w(23)=0.048048399081180627315979312780662;
      w(24)=0.033732303685955999377522748582756;
      w(25)=0.018896858024263465581352047676267;
      w(26)=0.0030769230769230769230769230769231;
   case(27)
      w(1)=0.0028490028490028490028490028490028;
      w(2)=0.017501974876065579019369734362719;
      w(3)=0.031262951735202384324784990064795;
      w(4)=0.044577657933061698744074742241386;
      w(5)=0.057265569680162731739087228253941;
      w(6)=0.069149342360043276280634735030127;
      w(7)=0.080062321970538458168238070483950;
      w(8)=0.089851365259290559972136301946063;
      w(9)=0.098379074585952763179308931269989;
      w(10)=0.10552574782125301121838198947990;
      w(11)=0.11119106525743703292826195934215;
      w(12)=0.11529550025465198281375333337473;
      w(13)=0.11778143658595615906649907944836;
      w(14)=0.11861397766276302708523980370575;
      w(15)=0.11778143658595615906649907944836;
      w(16)=0.11529550025465198281375333337473;
      w(17)=0.11119106525743703292826195934215;
      w(18)=0.10552574782125301121838198947990;
      w(19)=0.098379074585952763179308931269989;
      w(20)=0.089851365259290559972136301946063;
      w(21)=0.080062321970538458168238070483950;
      w(22)=0.069149342360043276280634735030127;
      w(23)=0.057265569680162731739087228253941;
      w(24)=0.044577657933061698744074742241386;
      w(25)=0.031262951735202384324784990064795;
      w(26)=0.017501974876065579019369734362719;
      w(27)=0.0028490028490028490028490028490028;
   case(28)
      w(1)=0.0026455026455026455026455026455026;
      w(2)=0.016255883957504218198990152897110;
      w(3)=0.029054220677979144706852779823582;
      w(4)=0.041466915243006721037985848778627;
      w(5)=0.053338077047327488400325513471500;
      w(6)=0.064513658080354538390245693462478;
      w(7)=0.074848123509707702780567939143162;
      w(8)=0.084206795121510325739491187797920;
      w(9)=0.092467685997712175722035394699386;
      w(10)=0.099523110412495672875945354750915;
      w(11)=0.10528109376105589448237755179491;
      w(12)=0.10966657379597662508858586379121;
      w(13)=0.11262238007723901254063642411366;
      w(14)=0.11410997967262783453331479283004;
      w(15)=0.11410997967262783453331479283004;
      w(16)=0.11262238007723901254063642411366;
      w(17)=0.10966657379597662508858586379121;
      w(18)=0.10528109376105589448237755179491;
      w(19)=0.099523110412495672875945354750915;
      w(20)=0.092467685997712175722035394699386;
      w(21)=0.084206795121510325739491187797920;
      w(22)=0.074848123509707702780567939143162;
      w(23)=0.064513658080354538390245693462478;
      w(24)=0.053338077047327488400325513471500;
      w(25)=0.041466915243006721037985848778627;
      w(26)=0.029054220677979144706852779823582;
      w(27)=0.016255883957504218198990152897110;
      w(28)=0.0026455026455026455026455026455026;
   case(29)
      w(1)=0.0024630541871921182266009852216749;
      w(2)=0.015138169859967596790857926606856;
      w(3)=0.027070806296824827540851425445026;
      w(4)=0.038668439979712979747832574267124;
      w(5)=0.049795809093237560922668341175478;
      w(6)=0.060318503828522735028064463517634;
      w(7)=0.070108938000597999185266858424014;
      w(8)=0.079048313027885601784170232175109;
      w(9)=0.087028133441135605637626155215633;
      w(10)=0.093951542114796312497836506987570;
      w(11)=0.099734501637149860653569177534068;
      w(12)=0.10430681646367653742270711219124;
      w(13)=0.10761298583356853937358685404661;
      w(14)=0.10961287784589042013213241075660;
      w(15)=0.11028221677968261011245795287074;
      w(16)=0.10961287784589042013213241075660;
      w(17)=0.10761298583356853937358685404661;
      w(18)=0.10430681646367653742270711219124;
      w(19)=0.099734501637149860653569177534068;
      w(20)=0.093951542114796312497836506987570;
      w(21)=0.087028133441135605637626155215633;
      w(22)=0.079048313027885601784170232175109;
      w(23)=0.070108938000597999185266858424014;
      w(24)=0.060318503828522735028064463517634;
      w(25)=0.049795809093237560922668341175478;
      w(26)=0.038668439979712979747832574267124;
      w(27)=0.027070806296824827540851425445026;
      w(28)=0.015138169859967596790857926606856;
      w(29)=0.0024630541871921182266009852216749;
   case(30)
      w(1)=0.0022988505747126436781609195402299;
      w(2)=0.014131799327905387640732168668208;
      w(3)=0.025283166740551402204268253850042;
      w(4)=0.036142094199408535314732683123360;
      w(5)=0.046590694533142927401880491446022;
      w(6)=0.056511197923080383302193710472303;
      w(7)=0.065791336397790054944101374593132;
      w(8)=0.074326003324718253834067642000831;
      w(9)=0.082018512833406914799721545527160;
      w(10)=0.088781712319765210167256434455921;
      w(11)=0.094538975193860891781132715364036;
      w(12)=0.099225071004299830657675964954309;
      w(13)=0.10278690530723498947028080833538;
      w(14)=0.10518412159645464985621298085914;
      w(15)=0.10638955872366792494758230680993;
      w(16)=0.10638955872366792494758230680993;
      w(17)=0.10518412159645464985621298085914;
      w(18)=0.10278690530723498947028080833538;
      w(19)=0.099225071004299830657675964954309;
      w(20)=0.094538975193860891781132715364036;
      w(21)=0.088781712319765210167256434455921;
      w(22)=0.082018512833406914799721545527160;
      w(23)=0.074326003324718253834067642000831;
      w(24)=0.065791336397790054944101374593132;
      w(25)=0.056511197923080383302193710472303;
      w(26)=0.046590694533142927401880491446022;
      w(27)=0.036142094199408535314732683123360;
      w(28)=0.025283166740551402204268253850042;
      w(29)=0.014131799327905387640732168668208;
      w(30)=0.0022988505747126436781609195402299;
   case(31)
      w(1)=0.0021505376344086021505376344086022;
      w(2)=0.013222471025464670302635629869835;
      w(3)=0.023666433230270316730045782709426;
      w(4)=0.033853940405224057629810846931873;
      w(5)=0.043681818160066912829747288083268;
      w(6)=0.053046465493448782774280860067640;
      w(7)=0.061848741290454623390748584002232;
      w(8)=0.069995377594100570450215353373140;
      w(9)=0.077400032341475618634283878417423;
      w(10)=0.083984220517529730665258108226376;
      w(11)=0.089678151045260822248142932836209;
      w(12)=0.094421468377857957432190489816534;
      w(13)=0.098163893013712757119741879401421;
      w(14)=0.10086575479865051339556259610285;
      w(15)=0.10249841359547039767590565499043;
      w(16)=0.10304456295320733314178496152549;
      w(17)=0.10249841359547039767590565499043;
      w(18)=0.10086575479865051339556259610285;
      w(19)=0.098163893013712757119741879401421;
      w(20)=0.094421468377857957432190489816534;
      w(21)=0.089678151045260822248142932836209;
      w(22)=0.083984220517529730665258108226376;
      w(23)=0.077400032341475618634283878417423;
      w(24)=0.069995377594100570450215353373140;
      w(25)=0.061848741290454623390748584002232;
      w(26)=0.053046465493448782774280860067640;
      w(27)=0.043681818160066912829747288083268;
      w(28)=0.033853940405224057629810846931873;
      w(29)=0.023666433230270316730045782709426;
      w(30)=0.013222471025464670302635629869835;
      w(31)=0.0021505376344086021505376344086022
   case(32)
      w(1)=0.0020161290322580645161290322580645;
      w(2)=0.012398106501373843788620349229117;
      w(3)=0.022199552889291964623832171161080;
      w(4)=0.031775135410915465781562278917906;
      w(5)=0.041034201586062723330403841719016;
      w(6)=0.049885271336221207011960153724443;
      w(7)=0.058240497248055869550798929919204;
      w(8)=0.066016877257154543932436331683494;
      w(9)=0.073137139602679032640370983569209;
      w(10)=0.079530525692106252292356728883049;
      w(11)=0.085133497949668230527527658506617;
      w(12)=0.089890372957357833072124789584664;
      w(13)=0.093753875546813813565908354145794;
      w(14)=0.096685608948002600560378147706395;
      w(15)=0.098656436540761777170651164243408;
      w(16)=0.099646771501276777634939084748541;
      w(17)=0.099646771501276777634939084748541;
      w(18)=0.098656436540761777170651164243408;
      w(19)=0.096685608948002600560378147706395;
      w(20)=0.093753875546813813565908354145794;
      w(21)=0.089890372957357833072124789584664;
      w(22)=0.085133497949668230527527658506617;
      w(23)=0.079530525692106252292356728883049;
      w(24)=0.073137139602679032640370983569209;
      w(25)=0.066016877257154543932436331683494;
      w(26)=0.058240497248055869550798929919204;
      w(27)=0.049885271336221207011960153724443;
      w(28)=0.041034201586062723330403841719016;
      w(29)=0.031775135410915465781562278917906;
      w(30)=0.022199552889291964623832171161080;
      w(31)=0.012398106501373843788620349229117;
      w(32)=0.0020161290322580645161290322580645;
   case(33)
      w(1)= 0.0018939393939393939393939393939394;
      w(2)=0.011648448392267734651222179870285;
      w(3)=0.020864609017603360095811664182613;
      w(4)=0.029881045916746477519971156454995;
      w(5)=0.038617814771813967563858875675217;
      w(6)=0.046993850461024170547973112041225;
      w(7)=0.054931059442626967951698841306657;
      w(8)=0.062355367852465305441060959167732;
      w(9)=0.069197469494016147559969796800326;
      w(10)=0.075393486923973828507129918402909;
      w(11)=0.080885572193455092179707696772196;
      w(12)=0.085622448531813132550191963745865;
      w(13)=0.089559889747077400661148462725561;
      w(14)=0.092661133442241463529935120285105;
      w(15)=0.094897224394591815824355414662604;
      w(16)=0.096247284972985461996554494277935;
      w(17)=0.096698710102716558960032808469672;
      w(18)=0.096247284972985461996554494277935;
      w(19)=0.094897224394591815824355414662604;
      w(20)=0.092661133442241463529935120285105;
      w(21)=0.089559889747077400661148462725561;
      w(22)=0.085622448531813132550191963745865;
      w(23)=0.080885572193455092179707696772196;
      w(24)=0.075393486923973828507129918402909;
      w(25)=0.069197469494016147559969796800326;
      w(26)=0.062355367852465305441060959167732;
      w(27)=0.054931059442626967951698841306657;
      w(28)=0.046993850461024170547973112041225;
      w(29)=0.038617814771813967563858875675217;
      w(30)=0.029881045916746477519971156454995;
      w(31)=0.020864609017603360095811664182613;
      w(32)=0.011648448392267734651222179870285;
      w(33)=0.0018939393939393939393939393939394;
   otherwise
      error('LOBATTO_WEIGHTS - Fatal error! Illegal value of n.');
end

return
end
