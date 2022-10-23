package riemann;

import java.math.BigDecimal;
import java.math.MathContext;
import java.math.RoundingMode;

/*
Author: O. Shanker.
oshanker At gmail dot com
http://sites.google.com/site/oshanker/

Calculate gram points for large t.

# Permission is hereby granted, free of charge, to any person obtaining 
# a copy of this software and associated documentation files (the
# "Software"), to deal in the Software without restriction, including
# without limitation the rights to use, copy, modify, merge, publish,
# distribute copies of the Software, and to
# permit persons to whom the Software is furnished to do so, subject to
# the following conditions:
# 
# The above copyright notice and this permission notice shall be
# included in all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
# MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
# IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
# CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
# TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
# SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/
public class Gram {

	static int[] primes = {2,3,5,7,11,13,17,19};
	public static MathContext mc = new MathContext(45, RoundingMode.HALF_EVEN);
	public static BigDecimal pi   = new BigDecimal("3.141592653589793238462643383279502884197169399");
	static BigDecimal e           = new BigDecimal("2.718281828459045235360287471352662497757247094");
	public static BigDecimal pie2 = new BigDecimal("17.07946844534713413092710173909314899006977707");
	public static BigDecimal bdTWO = BigDecimal.valueOf(2);
	public static BigDecimal pi8 = pi.divide(BigDecimal.valueOf(8), mc);
	static BigDecimal pi2 = pi.divide(bdTWO, mc);
	static BigDecimal pi4 = pi.divide(BigDecimal.valueOf(4), mc);
	public static BigDecimal pi_2 = pi.multiply(bdTWO, mc);	
    static BigDecimal iBD[] = new BigDecimal[20];
    static{
        for (int i = 0; i < 20; i++) {
            iBD[i] = BigDecimal.valueOf(i);
        }
    }
	
	/**
	 * High precision log values for the integers.
	 */
	public static BigDecimal[] logVals = {
			BigDecimal.ZERO,
		      new BigDecimal("0.693147180559945309417232121458176568075500134360255254120680009493393", mc),
		      new BigDecimal("1.098612288668109691395245236922525704647490557822749451734694333637494", mc),
		      new BigDecimal("1.386294361119890618834464242916353136151000268720510508241360018986787", mc),
		      new BigDecimal("1.609437912434100374600759333226187639525601354268517721912647891474178", mc),
		      new BigDecimal("1.791759469228055000812477358380702272722990692183004705855374343130887", mc),
		      new BigDecimal("1.945910149055313305105352743443179729637084729581861188459390149937579", mc),
		      new BigDecimal("2.079441541679835928251696364374529704226500403080765762362040028480180", mc),
		      new BigDecimal("2.197224577336219382790490473845051409294981115645498903469388667274988", mc),
		      new BigDecimal("2.302585092994045684017991454684364207601101488628772976033327900967572", mc),
		      new BigDecimal("2.397895272798370544061943577965129299821706853937417175218567709130573", mc),
		      new BigDecimal("2.484906649788000310229709479838878840798490826543259959976054352624281", mc),
		      new BigDecimal("2.564949357461536736053487441565318604805267944760207116419045510663464", mc),
		      new BigDecimal("2.639057329615258614522584864901356297712584863942116442580070159430973", mc),
		      new BigDecimal("2.708050201102210065996004570148713344173091912091267173647342225111673", mc),
		      new BigDecimal("2.772588722239781237668928485832706272302000537441021016482720037973574", mc),
		      new BigDecimal("2.833213344056216080249534617873126535588203012585744787297237737882292", mc),
		      new BigDecimal("2.890371757896164692207722595303227977370481250005754157590068676768382", mc),
		      new BigDecimal("2.944438979166440460009027431887853537237379261299128818537960236409292", mc),
		      new BigDecimal("2.995732273553990993435223576142540775676601622989028230154007910460966", mc),
		      new BigDecimal("3.044522437723422996500597980365705434284575287404610640194084483575074", mc),
		      new BigDecimal("3.091042453358315853479175699423305867897206988297672429339247718623967", mc),
		      new BigDecimal("3.135494215929149690806752831810196118442380314840435741998635377482993", mc),
		      new BigDecimal("3.178053830347945619646941601297055408873990960903515214096734362117675", mc),
		      new BigDecimal("3.218875824868200749201518666452375279051202708537035443825295782948357", mc),
		      new BigDecimal("3.258096538021482045470719563023495172880768079120462370539725520156858", mc),
		      new BigDecimal("3.295836866004329074185735710767577113942471673468248355204083000912482", mc),
		      new BigDecimal("3.332204510175203923939816986359532865788084998302371696700750168924367", mc),
		      new BigDecimal("3.367295829986474027183272032361911605494512913922744078921670351642780", mc),
		      new BigDecimal("3.401197381662155375413236691606889912248592046451522427768022234605066", mc),
		      new BigDecimal("3.433987204485146245929164324542357210449938930480591971756718072474981", mc),
		      new BigDecimal("3.465735902799726547086160607290882840377500671801276270603400047466968", mc),
		      new BigDecimal("3.496507561466480235457188814887655004469197411760166626953262042768067", mc),
		      new BigDecimal("3.526360524616161389666766739331303103663703146946000041417917747375686", mc),
		      new BigDecimal("3.555348061489413679706112076669367369162686083850378910372038041411758", mc),
		      new BigDecimal("3.583518938456110001624954716761404545445981384366009411710748686261775", mc),
		      new BigDecimal("3.610917912644224444368095671031447163900077587167636163644912681192989", mc),
		      new BigDecimal("3.637586159726385769426259553346030105312879395659384072658640245902686", mc),
		      new BigDecimal("3.663561646129646427448732678487844309452758502582956568153739844300958", mc),
		      new BigDecimal("3.688879454113936302852455697600717343752101757349283484274687919954359", mc),
		      new BigDecimal("3.713572066704307803866763373037407588376410469399301633619262910259978", mc),
		      new BigDecimal("3.737669618283368305917830101823882002360075421764865894314764493068467", mc),
		      new BigDecimal("3.761200115693562423472842513345847035559136184881555415191685264922859", mc),
		      new BigDecimal("3.784189633918261162896407820881482435972707122657927683459927728117360", mc),
		      new BigDecimal("3.806662489770319757391249807071239048820582469914016625382036558749167", mc),
		      new BigDecimal("3.828641396489095000223984953268372686517880449200690996119315386976386", mc),
		      new BigDecimal("3.850147601710058586820950669772173708896050502020224033200508346806818", mc),
		      new BigDecimal("3.871201010907890929064173722755231976949491095263770468217414371611068", mc),
		      new BigDecimal("3.891820298110626610210705486886359459274169459163722376918780299875159", mc),
		      new BigDecimal("3.912023005428146058618750787910551847126702842897290697945975792441751", mc),
		      new BigDecimal("3.931825632724325771644779854795652240235693570408494239031932071519786", mc),
		      new BigDecimal("3.951243718581427354887951684481671740956268213480717624660405529650251", mc),
		      new BigDecimal("3.970291913552121834144469139029057770359977752911217603048129470018004", mc),
		      new BigDecimal("3.988984046564274383602967832225753682017971807828503609324763010405876", mc),
		      new BigDecimal("4.007333185232470918662702911191316939347308208205934897131215600604752", mc),
		      new BigDecimal("4.025351690735149233357049107817709433863585132662626950821430178417760", mc),
		      new BigDecimal("4.043051267834550151404272668810379241884869819121878270272654570046786", mc),
		      new BigDecimal("4.060443010546419336600504153820088173570013048282999333042350361136174", mc),
		      new BigDecimal("4.077537443905719450616050373719697624063346789330454529512036697059200", mc),
		      new BigDecimal("4.094344562222100684830468813065066480324092180811777681888702244098460", mc),
		      new BigDecimal("4.110873864173311248751389103425614746315681743081261062937383646419439", mc),
		      new BigDecimal("4.127134385045091555346396446000533778525439064840847225877398081968375", mc),
		      new BigDecimal("4.143134726391532687895843217288231138932065845227360091928778817212568", mc),
		      new BigDecimal("4.158883083359671856503392728749059408453000806161531524724080056960361", mc),
		      new BigDecimal("4.174387269895637110654246774791506244330869299028724838331693402137643", mc),
		      new BigDecimal("4.189654742026425544874420936345831572544697546120421881073942052261461", mc),
		      new BigDecimal("4.204692619390966059670071996363722750566932903221895337137784130775268", mc),
		      new BigDecimal("4.219507705176106699083998860789479671739203281306255295538597756869079", mc),
		      new BigDecimal("4.234106504597259382201998068732721823089870872663185193733329711120487", mc),
		      new BigDecimal("4.248495242049358989123344198127543937238186218210634164492718050905152", mc),
		      new BigDecimal("4.262679877041315421329454532513034096759576526710566108121425802027351", mc),
		      new BigDecimal("4.276666119016055311042186838219581113521481518726264665831428695755169", mc),
		      new BigDecimal("4.290459441148391129092108857438542570904752844871597664595698857161789", mc),
		      new BigDecimal("4.304065093204169753785327792489623731975577721527891417765592690686383", mc),
		      new BigDecimal("4.317488113536310440596763903374900983698693266359784895559990116585852", mc),
		      new BigDecimal("4.330733340286331078843491674804206673388379530019639326779320255396079", mc),
		      new BigDecimal("4.343805421853683849167296321408309029458791583519278363677957859068153", mc),
		      new BigDecimal("4.356708826689591736865964799946020877528258636943211822274419853794352", mc),
		      new BigDecimal("4.369447852467021494172945541481410922173541224422609625412171117559806", mc),
		      new BigDecimal("4.382026634673881612269687819058893911827601891709538738395367929447753", mc),
		      new BigDecimal("4.394449154672438765580980947690102818589962231290997806938777334549977", mc),
		      new BigDecimal("4.406719247264253113283995494495584156451910603759556887739942919753372", mc),
		      new BigDecimal("4.418840607796597923475472223291370453029313056663236370187943462938578", mc),
		      new BigDecimal("4.43081679884331361533506222328205857043557555612512114843544450256186139991", mc),
		      new BigDecimal("4.442651256490316454850293951099314175113804366854262509209885629356471", mc),
		      new BigDecimal("4.454347296253507732890074634804023603634636319241810669312365274416252", mc),
		      new BigDecimal("4.465908118654583718578517269284437310142003471745493530656364685280275", mc),
		      new BigDecimal("4.477336814478206472313639942339659004048207257018182937580607737610754", mc),
		      new BigDecimal("4.488636369732139838317815540669849219404660387132959364106697577287953", mc),
		      new BigDecimal("4.499809670330265066808481928529415616896082604274271879502716568242561", mc),
		      new BigDecimal("4.510859506516850041158840185008498334442352674342068304878435660601044", mc),
		      new BigDecimal("4.521788577049040309641217074726549254593380583560946250239995396469780", mc),
		      new BigDecimal("4.532599493153255937324409561464882915097429488303341423491412406112475", mc),
		      new BigDecimal("4.543294782270003896238182791230350276971550636380479287321188356300211", mc),
		      new BigDecimal("4.553876891600540834609786765114041176762980615567646540450608127883471", mc),
		      new BigDecimal("4.564348191467836238481405844213408545024991229624025722338094381104462", mc),
		      new BigDecimal("4.574710978503382822116721621703961713808914902658781355976234368760177", mc),
		      new BigDecimal("4.584967478670571919627937608344536027349669593523977631039460309368553", mc),
		      new BigDecimal("4.595119850134589926852434051810180709116687969582916078687956376405562", mc),
		      new BigDecimal("4.605170185988091368035982909368728415202202977257545952066655801935145", mc),
	};
		

	public static void initLogVals(int N) {
		if (logVals.length < N) {
			BigDecimal[] temp = new BigDecimal[N];
			int oldlen = logVals.length;
			System.arraycopy(logVals, 0, temp, 0, oldlen);
			logVals = temp;
		}
	}
	
	static double logdbl(long arg) {
		return logVals[(int) (arg-1)].doubleValue();
	}
	
	public static BigDecimal log(long arg) {
		BigDecimal ret = null;
		if(arg<=logVals.length) {
		    ret = logVals[(int) (arg-1)];
			if(ret != null) {
				   return ret;
			}
		}
		for(int i = 0; i < primes.length; i++) {
		      if(arg%primes[i]==0) { 
		    	  ret = (logVals[primes[i]-1].add(log(arg/primes[i]), mc));
		  		if(arg<=logVals.length) {
		    	    logVals[(int) (arg-1)] = ret;
		  		}
			      return ret;
		      }
		}

		ret = (log(arg-1).add(
				log1sym(BigDecimal.ONE.divide(new BigDecimal(2*arg-1), mc ), mc), mc));
		if(arg<=logVals.length) {
		   logVals[(int) (arg-1)] = ret;
		}
	    return ret;
	}
	
	public static BigDecimal log(BigDecimal x, MathContext mc) {
		BigDecimal[] breakup = x.divideAndRemainder(BigDecimal.ONE);
		BigDecimal bigDecimalArg = breakup[0];
		BigDecimal delta = breakup[1];

		return log(bigDecimalArg.longValue()).add(
				log1sym(delta.divide(
						bigDecimalArg.add(x, mc), mc),mc),mc);
	}

	private static BigDecimal log1sym(BigDecimal x, MathContext mc) {
		BigDecimal term = x;
		BigDecimal mult = x.pow(2, mc);
		BigDecimal sum = x;
	    for (int i = 3; i < 20; i+= 2) {
	       term = term.multiply(mult, mc);
	       if(Math.abs(term.doubleValue()) < 1.0E-35){ 
	    	   break; 
	       }
	       //BigDecimal divisor = BigDecimal.valueOf(i);
	       sum = (sum.add(term.divide(iBD[i], mc),mc));
	    }
	    return sum.multiply(bdTWO, mc);
	}

	public static BigDecimal sqrt(BigDecimal x, MathContext mc, double prec) {
        double init = Math.sqrt(x.doubleValue());
        BigDecimal next = BigDecimal.valueOf(init);
        BigDecimal diff = null;
        while (true) {
            BigDecimal oldnext = next;
            next = x.divide(next, mc).add(next, mc).divide(Gram.bdTWO,mc);
            diff =  next.subtract(oldnext, mc);
            if(Math.abs(diff.doubleValue()) < prec){ 
                break; 
            }
        }
        return next;
	}
	
	public static BigDecimal theta(BigDecimal t, MathContext mc) {
		BigDecimal t2 = t.divide(bdTWO, mc);
		BigDecimal arg1 = sqrt(t2.divide(pi, mc), mc, 1.0E-18);
		initLogVals(arg1.intValue());
		BigDecimal term = BigDecimal.ONE.divide(t.multiply(BigDecimal.valueOf(48),mc), mc);
		BigDecimal theta = t.multiply(log(arg1,mc),mc).subtract(t2,mc).
				subtract(pi8 ,mc).add(term,mc);
		return theta;
	}	

	public static double thetaNormalized(BigDecimal t, MathContext mc){
		BigDecimal t2 = t.divide(Gram.bdTWO, mc);
		BigDecimal arg1 = t2.divide(Gram.pi, mc);
		BigDecimal sqrtArg1 = Gram.sqrt(arg1, Gram.mc, 1.0E-21);
		double theta = t.multiply(Gram.log(sqrtArg1, mc), mc).subtract(t2, mc)
				.subtract(Gram.pi8, mc).remainder(Gram.pi_2).doubleValue();
		return theta;
	}

	public static double gram(BigDecimal offset, double zero){
        BigDecimal tvalsi = offset.add(BigDecimal.valueOf(zero), Gram.mc);
        double diff = gramInterval(tvalsi);
        double heck1 = thetaNormalized(tvalsi, Gram.mc)  / Math.PI;
        if(heck1>=1){heck1--;}
        double gram0 = zero - (heck1*diff);
        return gram0;
	}

    public static long gramIndex(BigDecimal offset, double zero){
        BigDecimal tvalsi = offset.add(BigDecimal.valueOf(zero), Gram.mc);
        double heck1 = theta(tvalsi, Gram.mc).divide(pi, Gram.mc).doubleValue();
        return (long)(heck1+0.5);
    }
	
	public static void main(String[] args) throws Exception {
		test1();
	}

	private static void test1() {
//		BigDecimal offset = new BigDecimal("100000000000000000000000000");
//		BigDecimal zero = new BigDecimal("54.221610206411131994");
//		System.out.println(pi.multiply(e).multiply(bdTWO));
		BigDecimal offset = new BigDecimal("10000000000000000000000000000");
		BigDecimal zero = new BigDecimal("100.437512887104287873");
		BigDecimal zero1 = new BigDecimal("100.464843234223048518");
		BigDecimal tvalsi = offset.add(zero, Gram.mc);
		double diff = gramInterval(tvalsi);
		double heck1 = thetaNormalized(tvalsi, Gram.mc) *(diff)/Math.PI;
		
		double gram0 = zero.doubleValue() - (heck1);
        System.out.println(gram0 + ", cf riemann.ConjecturesTest.testX() " 
		+ (100.36777243017612) + " diff " + diff);
	}

    public static double gramInterval(BigDecimal tvalsi) {
        BigDecimal t2 = tvalsi.divide(bdTWO, mc);
		BigDecimal arg1 = sqrt(t2.divide(pi, mc), mc, 1.0E-15);
		BigDecimal fourthrootArg1 = Gram.sqrt(arg1, mc, 1.0E-15);
		Gram.initLogVals(fourthrootArg1.intValue());
		double diff = pi.divide(log(fourthrootArg1, mc).multiply(Gram.bdTWO), mc).doubleValue();
        return diff;
    }
	
	/*
awk '{print "new BigDecimal(\"" $2 ")"$3}' xx | sed -e 's/\(l..[1-9]\)/\1./' -e 's/,/",/'	 

*/

}
