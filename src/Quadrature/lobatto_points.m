function xi = lobatto_points(n)

xi = ones(1,n);

switch n
   case(2)
      xi(1)=-1.0000000000000000000000000000000;
      xi(2)=1.0000000000000000000000000000000;
   case(3)
      xi(1)=-1.0000000000000000000000000000000;
      xi(2)=0.0;
      xi(3)=1.0000000000000000000000000000000;
   case(4)
      xi(1)=-1.0000000000000000000000000000000;
      xi(2)=-0.44721359549995793928183473374626;
      xi(3)=0.44721359549995793928183473374626;
      xi(4)=1.0000000000000000000000000000000;
   case(5)
      xi(1)=-1.0000000000000000000000000000000;
      xi(2)=-0.65465367070797714379829245624686;
      xi(3)=0.0;
      xi(4)=0.65465367070797714379829245624686;
      xi(5)=1.0000000000000000000000000000000;
   case(6)
      xi(1)=-1.0000000000000000000000000000000;
      xi(2)=-0.76505532392946469285100297395934;
      xi(3)=-0.28523151648064509631415099404088;
      xi(4)=0.28523151648064509631415099404088;
      xi(5)=0.76505532392946469285100297395934;
      xi(6)=1.0000000000000000000000000000000;
   case(7)
      xi(1)=-1.0000000000000000000000000000000;
      xi(2)=-0.83022389627856692987203221396747;
      xi(3)=-0.46884879347071421380377188190877;
      xi(4)=0.0;
      xi(5)=0.46884879347071421380377188190877;
      xi(6)=0.83022389627856692987203221396747;
      xi(7)=1.0000000000000000000000000000000;
   case(8)
      xi(1)=-1.0000000000000000000000000000000;
      xi(2)=-0.87174014850960661533744576122066;
      xi(3)=-0.59170018143314230214451073139795;
      xi(4)=-0.20929921790247886876865726034535;
      xi(5)=0.20929921790247886876865726034535;
      xi(6)=0.59170018143314230214451073139795;
      xi(7)=0.87174014850960661533744576122066;
      xi(8)=1.0000000000000000000000000000000;
   case(9)
      xi(1)=-1.0000000000000000000000000000000;
      xi(2)=-0.89975799541146015731234524441834;
      xi(3)=-0.67718627951073775344588542709134;
      xi(4)=-0.36311746382617815871075206870866;
      xi(5)=0.0;
      xi(6)=0.36311746382617815871075206870866;
      xi(7)=0.67718627951073775344588542709134;
      xi(8)=0.89975799541146015731234524441834;
      xi(9)=1.0000000000000000000000000000000;
   case(10)
      xi(1)=-1.0000000000000000000000000000000;
      xi(2)=-0.91953390816645881382893266082234;
      xi(3)=-0.73877386510550507500310617485983;
      xi(4)=-0.47792494981044449566117509273126;
      xi(5)=-0.16527895766638702462621976595817;
      xi(6)=0.16527895766638702462621976595817;
      xi(7)=0.47792494981044449566117509273126;
      xi(8)=0.73877386510550507500310617485983;
      xi(9)=0.91953390816645881382893266082234;
      xi(10)=1.0000000000000000000000000000000;
   case(11)
      xi(1)=-1.0000000000000000000000000000000;
      xi(2)=-0.93400143040805913433227413609938;
      xi(3)=-0.78448347366314441862241781610846;
      xi(4)=-0.56523532699620500647096396947775;
      xi(5)=-0.29575813558693939143191151555906;
      xi(6)=0.0;
      xi(7)=0.29575813558693939143191151555906;
      xi(8)=0.56523532699620500647096396947775;
      xi(9)=0.78448347366314441862241781610846;
      xi(10)=0.93400143040805913433227413609938;
      xi(11)=1.0000000000000000000000000000000;
   case(12)
      xi(1)=-1.0000000000000000000000000000000;
      xi(2)=-0.94489927222288222340758013830322;
      xi(3)=-0.81927932164400667834864158171690;
      xi(4)=-0.63287615303186067766240485444366;
      xi(5)=-0.39953094096534893226434979156697;
      xi(6)=-0.13655293285492755486406185573969;
      xi(7)=0.13655293285492755486406185573969;
      xi(8)=0.39953094096534893226434979156697;
      xi(9)=0.63287615303186067766240485444366;
      xi(10)=0.81927932164400667834864158171690;
      xi(11)=0.94489927222288222340758013830322;
      xi(12)=1.0000000000000000000000000000000;
   case(13)
      xi(1)=-1.0000000000000000000000000000000;
      xi(2)=-0.95330984664216391189690546475545;
      xi(3)=-0.84634756465187231686592560709875;
      xi(4)=-0.68618846908175742607275903956636;
      xi(5)=-0.48290982109133620174693723363693;
      xi(6)=-0.24928693010623999256867370037423;
      xi(7)=0.0;
      xi(8)=0.24928693010623999256867370037423;
      xi(9)=0.48290982109133620174693723363693;
      xi(10)=0.68618846908175742607275903956636;
      xi(11)=0.84634756465187231686592560709875;
      xi(12)=0.95330984664216391189690546475545;
      xi(13)=1.0000000000000000000000000000000;
   case(14)
      xi(1)=-1.0000000000000000000000000000000;
      xi(2)=-0.95993504526726090135510016201542;
      xi(3)=-0.86780105383034725100022020290826;
      xi(4)=-0.72886859909132614058467240052088;
      xi(5)=-0.55063940292864705531662270585908;
      xi(6)=-0.34272401334271284504390340364167;
      xi(7)=-0.11633186888370386765877670973616;
      xi(8)=0.11633186888370386765877670973616;
      xi(9)=0.34272401334271284504390340364167;
      xi(10)=0.55063940292864705531662270585908;
      xi(11)=0.72886859909132614058467240052088;
      xi(12)=0.86780105383034725100022020290826;
      xi(13)=0.95993504526726090135510016201542;
      xi(14)=1.0000000000000000000000000000000;
   case(15)
      xi(1)=-1.0000000000000000000000000000000;
      xi(2)=-0.96524592650383857279585139206960;
      xi(3)=-0.88508204422297629882540163148223;
      xi(4)=-0.76351968995181520070411847597629;
      xi(5)=-0.60625320546984571112352993863673;
      xi(6)=-0.42063805471367248092189693873858;
      xi(7)=-0.21535395536379423822567944627292;
      xi(8)=0.0;
      xi(9)=0.21535395536379423822567944627292;
      xi(10)=0.42063805471367248092189693873858;
      xi(11)=0.60625320546984571112352993863673;
      xi(12)=0.76351968995181520070411847597629;
      xi(13)=0.88508204422297629882540163148223;
      xi(14)=0.96524592650383857279585139206960;
      xi(15)=1.0000000000000000000000000000000;
   case(16)
      xi(1)=-1.0000000000000000000000000000000;
      xi(2)=-0.96956804627021793295224273836746;
      xi(3)=-0.89920053309347209299462826151985;
      xi(4)=-0.79200829186181506393108827096315;
      xi(5)=-0.65238870288249308946788321964058;
      xi(6)=-0.48605942188713761178189078584687;
      xi(7)=-0.29983046890076320809835345472230;
      xi(8)=-0.10132627352194944784303300504592;
      xi(9)=0.10132627352194944784303300504592;
      xi(10)=0.29983046890076320809835345472230;
      xi(11)=0.48605942188713761178189078584687;
      xi(12)=0.65238870288249308946788321964058;
      xi(13)=0.79200829186181506393108827096315;
      xi(14)=0.89920053309347209299462826151985;
      xi(15)=0.96956804627021793295224273836746;
      xi(16)=1.0000000000000000000000000000000;
   case(17)
      xi(1)=-1.0000000000000000000000000000000;
      xi(2)=-0.97313217663141831415697950187372;
      xi(3)=-0.91087999591557359562380250639773;
      xi(4)=-0.81569625122177030710675055323753;
      xi(5)=-0.69102898062768470539491935737245;
      xi(6)=-0.54138539933010153912373340750406;
      xi(7)=-0.37217443356547704190723468073526;
      xi(8)=-0.18951197351831738830426301475311;
      xi(9)=0.0;
      xi(10)=0.18951197351831738830426301475311;
      xi(11)=0.37217443356547704190723468073526;
      xi(12)=0.54138539933010153912373340750406;
      xi(13)=0.69102898062768470539491935737245;
      xi(14)=0.81569625122177030710675055323753;
      xi(15)=0.91087999591557359562380250639773;
      xi(16)=0.97313217663141831415697950187372;
      xi(17)=1.0000000000000000000000000000000;
   case(18)
      xi(1)=-1.0000000000000000000000000000000;
      xi(2)=-0.97610555741219854286451892434170;
      xi(3)=-0.92064918534753387383785462543128;
      xi(4)=-0.83559353521809021371364636232794;
      xi(5)=-0.72367932928324268130621036530207;
      xi(6)=-0.58850483431866176117353589319356;
      xi(7)=-0.43441503691212397534228713674067;
      xi(8)=-0.26636265287828098416766533202560;
      xi(9)=-0.089749093484652111022645010088562;
      xi(10)=0.089749093484652111022645010088562;
      xi(11)=0.26636265287828098416766533202560;
      xi(12)=0.43441503691212397534228713674067;
      xi(13)=0.58850483431866176117353589319356;
      xi(14)=0.72367932928324268130621036530207;
      xi(15)=0.83559353521809021371364636232794;
      xi(16)=0.92064918534753387383785462543128;
      xi(17)=0.97610555741219854286451892434170;
      xi(18)=1.0000000000000000000000000000000;
   case(19)
      xi(1)=-1.0000000000000000000000000000000;
      xi(2)=-0.97861176622208009515263406311022;
      xi(3)=-0.92890152815258624371794025879655;
      xi(4)=-0.85246057779664609308595597004106;
      xi(5)=-0.75149420255261301416363748963394;
      xi(6)=-0.62890813726522049776683230622873;
      xi(7)=-0.48822928568071350277790963762492;
      xi(8)=-0.33350484782449861029850010384493;
      xi(9)=-0.16918602340928157137515415344488;
      xi(10)=0.0;
      xi(11)=0.16918602340928157137515415344488;
      xi(12)=0.33350484782449861029850010384493;
      xi(13)=0.48822928568071350277790963762492;
      xi(14)=0.62890813726522049776683230622873;
      xi(15)=0.75149420255261301416363748963394;
      xi(16)=0.85246057779664609308595597004106;
      xi(17)=0.92890152815258624371794025879655;
      xi(18)=0.97861176622208009515263406311022;
      xi(19)=1.0000000000000000000000000000000;
   case(20)
      xi(1)=-1.0000000000000000000000000000000;
      xi(2)=-0.98074370489391417192544643858423;
      xi(3)=-0.93593449881266543571618158493063;
      xi(4)=-0.86687797808995014130984721461629;
      xi(5)=-0.77536826095205587041431752759469;
      xi(6)=-0.66377640229031128984640332297116;
      xi(7)=-0.53499286403188626164813596182898;
      xi(8)=-0.39235318371390929938647470381582;
      xi(9)=-0.23955170592298649518240135692709;
      xi(10)=-0.080545937238821837975944518159554;
      xi(11)=0.080545937238821837975944518159554;
      xi(12)=0.23955170592298649518240135692709;
      xi(13)=0.39235318371390929938647470381582;
      xi(14)=0.53499286403188626164813596182898;
      xi(15)=0.66377640229031128984640332297116;
      xi(16)=0.77536826095205587041431752759469;
      xi(17)=0.86687797808995014130984721461629;
      xi(18)=0.93593449881266543571618158493063;
      xi(19)=0.98074370489391417192544643858423;
      xi(20)=1.0000000000000000000000000000000;
   case(21)
      xi(1)=-1.0000000000000000000000000000000;
      xi(2)=-0.98257229660454802823448127655541;
      xi(3)=-0.94197629695974553429610265066144;
      xi(4)=-0.87929475532359046445115359630494;
      xi(5)=-0.79600192607771240474431258966036;
      xi(6)=-0.69405102606222323262731639319467;
      xi(7)=-0.57583196026183068692702187033809;
      xi(8)=-0.44411578327900210119451634960735;
      xi(9)=-0.30198985650876488727535186785875;
      xi(10)=-0.15278551580218546600635832848567;
      xi(11)=0.0;
      xi(12)=0.15278551580218546600635832848567;
      xi(13)=0.30198985650876488727535186785875;
      xi(14)=0.44411578327900210119451634960735;
      xi(15)=0.57583196026183068692702187033809;
      xi(16)=0.69405102606222323262731639319467;
      xi(17)=0.79600192607771240474431258966036;
      xi(18)=0.87929475532359046445115359630494;
      xi(19)=0.94197629695974553429610265066144;
      xi(20)=0.98257229660454802823448127655541;
      xi(21)=1.0000000000000000000000000000000;
   case(22)
      xi(1)=-1.0000000000000000000000000000000;
      xi(2)=-0.98415243845764617655228962221207;
      xi(3)=-0.94720428399922868052421376661573;
      xi(4)=-0.89006229019090447052965782577909;
      xi(5)=-0.81394892761192113604544184805614;
      xi(6)=-0.72048723996120215811988189639847;
      xi(7)=-0.61166943828425897122621160586993;
      xi(8)=-0.48981487518990234980875123568327;
      xi(9)=-0.35752071013891953806095728024018;
      xi(10)=-0.21760658515928504178795509346539;
      xi(11)=-0.073054540010898334761088790464107;
      xi(12)=0.073054540010898334761088790464107;
      xi(13)=0.21760658515928504178795509346539;
      xi(14)=0.35752071013891953806095728024018;
      xi(15)=0.48981487518990234980875123568327;
      xi(16)=0.61166943828425897122621160586993;
      xi(17)=0.72048723996120215811988189639847;
      xi(18)=0.81394892761192113604544184805614;
      xi(19)=0.89006229019090447052965782577909;
      xi(20)=0.94720428399922868052421376661573;
      xi(21)=0.98415243845764617655228962221207;
      xi(22)=1.0000000000000000000000000000000;
   case(23)
      xi(1)=-1.0000000000000000000000000000000;
      xi(2)=-0.98552715587873257808146276673810;
      xi(3)=-0.95175795571071020413563967985143;
      xi(4)=-0.89945855804034501095016032034737;
      xi(5)=-0.82965109665128588622320061929000;
      xi(6)=-0.74369504117206068394516354306700;
      xi(7)=-0.64326364446013620847614553360277;
      xi(8)=-0.53031177113684416813011532015230;
      xi(9)=-0.40703793791447482919595048821510;
      xi(10)=-0.27584154894579306710687763267914;
      xi(11)=-0.13927620404066839859186261298277;
      xi(12)=0.0;
      xi(13)=0.13927620404066839859186261298277;
      xi(14)=0.27584154894579306710687763267914;
      xi(15)=0.40703793791447482919595048821510;
      xi(16)=0.53031177113684416813011532015230;
      xi(17)=0.64326364446013620847614553360277;
      xi(18)=0.74369504117206068394516354306700;
      xi(19)=0.82965109665128588622320061929000;
      xi(20)=0.89945855804034501095016032034737;
      xi(21)=0.95175795571071020413563967985143;
      xi(22)=0.98552715587873257808146276673810;
      xi(23)=1.0000000000000000000000000000000;
   case(24)
      xi(1)=-1.0000000000000000000000000000000;
      xi(2)=-0.98673055350516088355308673815447;
      xi(3)=-0.95574822092988635802697713055064;
      xi(4)=-0.90770567511350652199515299646621;
      xi(5)=-0.84346407015487204062330503742334;
      xi(6)=-0.76417048242049330778737528095229;
      xi(7)=-0.67124010526412869983566485818701;
      xi(8)=-0.56633135797929531218940954454228;
      xi(9)=-0.45131637321432261824821849156962;
      xi(10)=-0.32824761337551091203338917935961;
      xi(11)=-0.19932125339083266723657253912499;
      xi(12)=-0.066837993737228578113641808391677;
      xi(13)=0.066837993737228578113641808391677;
      xi(14)=0.19932125339083266723657253912499;
      xi(15)=0.32824761337551091203338917935961;
      xi(16)=0.45131637321432261824821849156962;
      xi(17)=0.56633135797929531218940954454228;
      xi(18)=0.67124010526412869983566485818701;
      xi(19)=0.76417048242049330778737528095229;
      xi(20)=0.84346407015487204062330503742334;
      xi(21)=0.90770567511350652199515299646621;
      xi(22)=0.95574822092988635802697713055064;
      xi(23)=0.98673055350516088355308673815447;
      xi(24)=1.0000000000000000000000000000000;
   case(25)
      xi(1)=-1.0000000000000000000000000000000;
      xi(2)=-0.98778994493149370927180407087316;
      xi(3)=-0.95926413825253447885984517883952;
      xi(4)=-0.91498277073462257832314993373685;
      xi(5)=-0.85567646583531657752381326017175;
      xi(6)=-0.78231965924071678039918795222875;
      xi(7)=-0.69611704881513436676036543378938;
      xi(8)=-0.59848414727999326809762099718107;
      xi(9)=-0.49102411481887838261895992256880;
      xi(10)=-0.37550145785922723322871461226072;
      xi(11)=-0.25381306416887658017988688168728;
      xi(12)=-0.12795705948310697270898462509465;
      xi(13)=0.0;
      xi(14)=0.12795705948310697270898462509465;
      xi(15)=0.25381306416887658017988688168728;
      xi(16)=0.37550145785922723322871461226072;
      xi(17)=0.49102411481887838261895992256880;
      xi(18)=0.59848414727999326809762099718107;
      xi(19)=0.69611704881513436676036543378938;
      xi(20)=0.78231965924071678039918795222875;
      xi(21)=0.85567646583531657752381326017175;
      xi(22)=0.91498277073462257832314993373685;
      xi(23)=0.95926413825253447885984517883952;
      xi(24)=0.98778994493149370927180407087316;
      xi(25)=1.0000000000000000000000000000000;
   case(26)
      xi(1)=-1.0000000000000000000000000000000;
      xi(2)=-0.98872741231147565442924669037130;
      xi(3)=-0.96237787476771732975066267128213;
      xi(4)=-0.92143554681755734112513052174869;
      xi(5)=-0.86652432395912357097445149780977;
      xi(6)=-0.79847718310743743903678763857034;
      xi(7)=-0.71832581636266508052540538790393;
      xi(8)=-0.62728529949231688780051890716122;
      xi(9)=-0.52673574202987854529824339438099;
      xi(10)=-0.41820138706624678556368388233385;
      xi(11)=-0.30332751285925272077404103790211;
      xi(12)=-0.18385549527005490132120567188576;
      xi(13)=-0.061596411781919728205487274781751;
      xi(14)=0.061596411781919728205487274781751;
      xi(15)=0.18385549527005490132120567188576;
      xi(16)=0.30332751285925272077404103790211;
      xi(17)=0.41820138706624678556368388233385;
      xi(18)=0.52673574202987854529824339438099;
      xi(19)=0.62728529949231688780051890716122;
      xi(20)=0.71832581636266508052540538790393;
      xi(21)=0.79847718310743743903678763857034;
      xi(22)=0.86652432395912357097445149780977;
      xi(23)=0.92143554681755734112513052174869;
      xi(24)=0.96237787476771732975066267128213;
      xi(25)=0.98872741231147565442924669037130;
      xi(26)=1.0000000000000000000000000000000;
   case(27)
      xi(1)=-1.0000000000000000000000000000000;
      xi(2)=-0.98956096372855062114887464005452;
      xi(3)=-0.96514840245081891926426680312113;
      xi(4)=-0.92718345872511577885974676407272;
      xi(5)=-0.87620208621452247866654304919414;
      xi(6)=-0.81292048689581228825809186475854;
      xi(7)=-0.73822714984645991083942722574841;
      xi(8)=-0.65317066369680951517870487518246;
      xi(9)=-0.55894506094256118006266807486812;
      xi(10)=-0.45687307561408240367861265083924;
      xi(11)=-0.34838758198902867043708657858717;
      xi(12)=-0.23501148310291813333680158051333;
      xi(13)=-0.11833633389852104843870826827532;
      xi(14)=0.0;
      xi(15)=0.11833633389852104843870826827532;
      xi(16)=0.23501148310291813333680158051333;
      xi(17)=0.34838758198902867043708657858717;
      xi(18)=0.45687307561408240367861265083924;
      xi(19)=0.55894506094256118006266807486812;
      xi(20)=0.65317066369680951517870487518246;
      xi(21)=0.73822714984645991083942722574841;
      xi(22)=0.81292048689581228825809186475854;
      xi(23)=0.87620208621452247866654304919414;
      xi(24)=0.92718345872511577885974676407272;
      xi(25)=0.96514840245081891926426680312113;
      xi(26)=0.98956096372855062114887464005452;
      xi(27)=1.0000000000000000000000000000000;
   case(28)
      xi(1)=-1.0000000000000000000000000000000;
      xi(2)=-0.99030540261845412551709291186460;
      xi(3)=-0.96762428585713130434377087914347;
      xi(4)=-0.93232516712155852492614241860899;
      xi(5)=-0.88487101721130283841230633825725;
      xi(6)=-0.82588097005633819604451442619660;
      xi(7)=-0.75612419400556976515209572526619;
      xi(8)=-0.67651012892957331771938124507339;
      xi(9)=-0.58807668983717560647831207976171;
      xi(10)=-0.49197675393157938042961668294448;
      xi(11)=-0.38946313757636280793588415520411;
      xi(12)=-0.28187226662160237195415250890933;
      xi(13)=-0.17060675530800436088270267774877;
      xi(14)=-0.057117121693512897626543651573747;
      xi(15)=0.057117121693512897626543651573747;
      xi(16)=0.17060675530800436088270267774877;
      xi(17)=0.28187226662160237195415250890933;
      xi(18)=0.38946313757636280793588415520411;
      xi(19)=0.49197675393157938042961668294448;
      xi(20)=0.58807668983717560647831207976171;
      xi(21)=0.67651012892957331771938124507339;
      xi(22)=0.75612419400556976515209572526619;
      xi(23)=0.82588097005633819604451442619660;
      xi(24)=0.88487101721130283841230633825725;
      xi(25)=0.93232516712155852492614241860899;
      xi(26)=0.96762428585713130434377087914347;
      xi(27)=0.99030540261845412551709291186460;
      xi(28)=1.0000000000000000000000000000000;
   case(29)
      xi(1)=-1.0000000000000000000000000000000;
      xi(2)=-0.99097298826856981627874511329957;
      xi(3)=-0.96984580728793629371240138788782;
      xi(4)=-0.93694271852098251154203684828336;
      xi(5)=-0.89266571997608827876422991405309;
      xi(6)=-0.83755273628178656524207331618682;
      xi(7)=-0.77227289720646874135417063413645;
      xi(8)=-0.69761866135636798839397442857709;
      xi(9)=-0.61449625220343294997090361742436;
      xi(10)=-0.52391467437196908357881657312361;
      xi(11)=-0.42697347171349438763757641096836;
      xi(12)=-0.32484938284191100626001402307884;
      xi(13)=-0.21878205828426101980905185232499;
      xi(14)=-0.11005901339559210970778760225568;
      xi(15)=0.0;
      xi(16)=0.11005901339559210970778760225568;
      xi(17)=0.21878205828426101980905185232499;
      xi(18)=0.32484938284191100626001402307884;
      xi(19)=0.42697347171349438763757641096836;
      xi(20)=0.52391467437196908357881657312361;
      xi(21)=0.61449625220343294997090361742436;
      xi(22)=0.69761866135636798839397442857709;
      xi(23)=0.77227289720646874135417063413645;
      xi(24)=0.83755273628178656524207331618682;
      xi(25)=0.89266571997608827876422991405309;
      xi(26)=0.93694271852098251154203684828336;
      xi(27)=0.96984580728793629371240138788782;
      xi(28)=0.99097298826856981627874511329957;
      xi(29)=1.0000000000000000000000000000000;
   case(30)
      xi(1)=-1.0000000000000000000000000000000;
      xi(2)=-0.99157394284050029333881822590865;
      xi(3)=-0.97184660316626924167659295220919;
      xi(4)=-0.94110478095105708230718954941552;
      xi(5)=-0.89969921819927685955334282314945;
      xi(6)=-0.84809948718019810955142383218307;
      xi(7)=-0.78689035723754708044995952059240;
      xi(8)=-0.71676539863708513163379859612990;
      xi(9)=-0.63851917580755840737109348270974;
      xi(10)=-0.55303826009505285238455023004066;
      xi(11)=-0.46129119016824068522655722110715;
      xi(12)=-0.36431750042244899775598510555179;
      xi(13)=-0.26321594371957379126709946900007;
      xi(14)=-0.15913204262585046782494857709885;
      xi(15)=-0.053245110485486669363007834451354;
      xi(16)=0.053245110485486669363007834451354;
      xi(17)=0.15913204262585046782494857709885;
      xi(18)=0.26321594371957379126709946900007;
      xi(19)=0.36431750042244899775598510555179;
      xi(20)=0.46129119016824068522655722110715;
      xi(21)=0.55303826009505285238455023004066;
      xi(22)=0.63851917580755840737109348270974;
      xi(23)=0.71676539863708513163379859612990;
      xi(24)=0.78689035723754708044995952059240;
      xi(25)=0.84809948718019810955142383218307;
      xi(26)=0.89969921819927685955334282314945;
      xi(27)=0.94110478095105708230718954941552;
      xi(28)=0.97184660316626924167659295220919;
      xi(29)=0.99157394284050029333881822590865;
      xi(30)=1.0000000000000000000000000000000;
   case(31)
      xi(1)=-1.0000000000000000000000000000000;
      xi(2)=-0.99211684434648108858932461564871;
      xi(3)=-0.97365493581573647482606954847523;
      xi(4)=-0.94486917020803922980640553764986;
      xi(5)=-0.90606695144126981314451152324746;
      xi(6)=-0.85765999529745561036700538239647;
      xi(7)=-0.80016154319246298927233133398318;
      xi(8)=-0.73418113630907524213414083483764;
      xi(9)=-0.66041820261152552333704127459677;
      xi(10)=-0.57965465720800191300329686219773;
      xi(11)=-0.49274661909883240570380023307021;
      xi(12)=-0.40061533828056180630266441939609;
      xi(13)=-0.30423743127287262095945514762626;
      xi(14)=-0.20463452924752450114109358572771;
      xi(15)=-0.10286244876068227708911472330805;
      xi(16)=0.0;
      xi(17)=0.10286244876068227708911472330805;
      xi(18)=0.20463452924752450114109358572771;
      xi(19)=0.30423743127287262095945514762626;
      xi(20)=0.40061533828056180630266441939609;
      xi(21)=0.49274661909883240570380023307021;
      xi(22)=0.57965465720800191300329686219773;
      xi(23)=0.66041820261152552333704127459677;
      xi(24)=0.73418113630907524213414083483764;
      xi(25)=0.80016154319246298927233133398318;
      xi(26)=0.85765999529745561036700538239647;
      xi(27)=0.90606695144126981314451152324746;
      xi(28)=0.94486917020803922980640553764986;
      xi(29)=0.97365493581573647482606954847523;
      xi(30)=0.99211684434648108858932461564871;
      xi(31)=1.0000000000000000000000000000000;
   case(32)
      xi(1)=-1.0000000000000000000000000000000;
      xi(2)=-0.99260893397276135937154131527745;
      xi(3)=-0.97529469048270922806245595038815;
      xi(4)=-0.94828483841723237808332824158843;
      xi(5)=-0.91184993906373190407453145237407;
      xi(6)=-0.86635247601267551983081476942454;
      xi(7)=-0.81224473177744234454682178807258;
      xi(8)=-0.75006449393667479771727596467830;
      xi(9)=-0.68042975561555081594239575769885;
      xi(10)=-0.60403258714842112613717291261189;
      xi(11)=-0.52163226288156529060659257339202;
      xi(12)=-0.43404771720184693960331935700639;
      xi(13)=-0.34214940653888148625381222069872;
      xi(14)=-0.24685065885020530441624143115887;
      xi(15)=-0.14909859681364749491438187617004;
      xi(16)=-0.049864725046593252306297744728426;
      xi(17)=0.049864725046593252306297744728426;
      xi(18)=0.14909859681364749491438187617004;
      xi(19)=0.24685065885020530441624143115887;
      xi(20)=0.34214940653888148625381222069872;
      xi(21)=0.43404771720184693960331935700639;
      xi(22)=0.52163226288156529060659257339202;
      xi(23)=0.60403258714842112613717291261189;
      xi(24)=0.68042975561555081594239575769885;
      xi(25)=0.75006449393667479771727596467830;
      xi(26)=0.81224473177744234454682178807258;
      xi(27)=0.86635247601267551983081476942454;
      xi(28)=0.91184993906373190407453145237407; 
      xi(29)=0.94828483841723237808332824158843; 
      xi(30)=0.97529469048270922806245595038815; 
      xi(31)=0.99260893397276135937154131527745; 
      xi(32)=1.0000000000000000000000000000000; case(33) 
      xi(1)=-1.0000000000000000000000000000000; 
      xi(2)=-0.99305635843365834366733608675624; 
      xi(3)=-0.97678616331690630148562760838132; 
      xi(4)=-0.95139345139699574337552441691356; 
      xi(5)=-0.91711730345094124082567906975682; 
      xi(6)=-0.87427810075056222064656062933034; 
      xi(7)=-0.82327592300406746955907131326671; 
      xi(8)=-0.76458700179352862800665848042309; 
      xi(9)=-0.69875931661816259569967724754442;
      xi(10)=-0.62640749128126825726383134241071;
      xi(11)=-0.54820705991911162311040471176063;
      xi(12)=-0.46488816163210675597300372098526;
      xi(13)=-0.37722872425339363504193466071800;
      xi(14)=-0.2860472014876740416053445029630;
      xi(15)=-0.19219493146747722574130699544116;
      xi(16)=-0.096548188176107006316957788091625;
      xi(17)=0.0;
      xi(18)=0.096548188176107006316957788091625;
      xi(19)=0.19219493146747722574130699544116;
      xi(20)=0.28604720148767404160534450296304;
      xi(21)=0.37722872425339363504193466071800;
      xi(22)=0.46488816163210675597300372098526;
      xi(23)=0.54820705991911162311040471176063;
      xi(24)=0.62640749128126825726383134241071;
      xi(25)=0.69875931661816259569967724754442;
      xi(26)=0.76458700179352862800665848042309;
      xi(27)=0.82327592300406746955907131326671;
      xi(28)=0.87427810075056222064656062933034;
      xi(29)=0.91711730345094124082567906975682;
      xi(30)=0.95139345139699574337552441691356;
      xi(31)=0.97678616331690630148562760838132;
      xi(32)=0.99305635843365834366733608675624;
      xi(33)=1.0000000000000000000000000000000;
   otherwise
      error('LOBATTO_POINTS - Fatal error! Illegal value of n.');
end

return
end
