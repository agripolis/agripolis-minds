/*************************************************************************
* This file is part of AgriPoliS-MINDS
*
* AgriPoliS: An Agricultural Policy Simulator
*
* Copyright (c) 2023 Alfons Balmann, Kathrin Happe, Konrad Kellermann et al.
* (cf. AUTHORS.md) at Leibniz Institute of Agricultural Development in 
* Transition Economies
*
* SPDX-License-Identifier: MIT
**************************************************************************/

// RegSurrogate.cpp
//---------------------------------------------------------------------------
#include <fstream>
#include <stdio.h>
#include <time.h>
#include <vector>
#include <sstream>
#include <iterator>
#include <istream>
#include <iomanip>
#include <filesystem>

#include "RegSurrogate.h"
#include "RegManager.h"
#include "RegPlot.h"
#include "surrogate.h"

#include "textinput.h"
namespace fs = std::filesystem;

enum objtype { PROD, INVEST, OBJVALUE, UNKNOWN };
static map <string, int> marketId;
static map < string, int> investId;
static pair<objtype, int> result_type(string aname) {
    objtype ot;
    int src;
    if (marketId.find(aname) != marketId.end()) {
        src = marketId[aname];
        ot = PROD;
    }
    else if (investId.find(aname) != investId.end()) {
        src = investId[aname];
        ot = INVEST;
    }
    else if (aname == "Objective_Value") {
        src = -1;
        ot = OBJVALUE;
    }
    else {
        src = -1;
        ot = UNKNOWN;
    }
    return pair<objtype, int>(ot, src);
}


static string rtrim(string s, char c) {
	int n = s.size();
	string rs,res;
	rs.resize(n);
	copy(s.rbegin(), s.rend(), rs.begin());
	size_t pos=0;
	while (rs[pos] == '\\') ++pos;
	
	res.resize(n - pos);
	for (unsigned i = 0; i < n-pos; ++i) {
		res[i] = rs[n - 1- i];
	}
	return res;
}

RegSurrogate::RegSurrogate() {
    obj_backup=NULL;
    ofstream out;
    flat_copy=false;
}

RegSurrogate::~RegSurrogate() {
    if (!flat_copy) {
        for (unsigned i=0;i<invest_links.size();i++)
            if (invest_links[i]) delete invest_links[i];
        for (unsigned i=0;i<market_links.size();i++)
            if (market_links[i])delete market_links[i];
        for (unsigned i=0;i<reference_links.size();i++)
            if (reference_links[i])delete reference_links[i];
        for (unsigned i=0;i<number_links.size();i++)
            if (number_links[i])delete number_links[i];
        for (unsigned i=0;i<land_links.size();i++)
            if (land_links[i]) delete land_links[i];

		for (unsigned i=0;i<yield_links.size();i++)
            if (yield_links[i]) delete yield_links[i];
    }
	yield_links.clear();

    invest_links.clear();
    market_links.clear();
    reference_links.clear();
    number_links.clear();
    land_links.clear();
    incomepay_links.clear();
    if (obj_backup) delete obj_backup;
}

void RegSurrogate::debug(string filename, bool before=false) {
    ofstream out;
    out.open(filename.c_str(),ios::trunc);
    out << "Inputs: " << endl;
    for (int i = 0; i < surrogateIO.inputlinks.size(); ++i) {
        out << surrogateIO.inputlinks[i].name << "\t" << inputs[i] << endl;
    }
    if (!before) {
        out << "\nOutputs: " << endl;
        for (int i = 0; i < surrogateIO.output_names.size(); ++i) {
            out << surrogateIO.output_names[i] << "\t" << x[i] << endl;
        }
    }
    out.close();
}

void RegSurrogate::debugPredict(string filename, vector<float> inp, vector<float> outp) {
    ofstream out;
    out.open((filename+"_in.csv").c_str(), ios::trunc);
    //out << "Inputs: " << endl;
    for (int i = 0; i < surrogateIO.inputlinks.size(); ++i) {
        out << surrogateIO.inputlinks[i].name << ";" << inp[i] << endl;
    }
    out.close();

    out.open((filename + "_out.csv").c_str(), ios::trunc);
    //out << "Inputs: " << endl;
    for (int i = 0; i < surrogateIO.output_names.size(); ++i) {
            out << surrogateIO.output_names[i] << ";" << outp[i] << endl;
    }
    out.close();
}

static string remove_price_suffix(string str) {
    string res=str;
    unsigned int len=res.length();
    if (len>5 && res.substr(len-5,5).compare("Price") == 0) {
        res = res.substr(0, len-5);
    }
    return res;
}

void RegSurrogate::debugCSV(vector<float> inps, string filename, bool before = true) {
    ofstream out;
    out.open((filename+"_in.csv").c_str(), ios::trunc);
    for (int i = 0; i < surrogateIO.inputlinks.size(); ++i) {
        out << surrogateIO.inputlinks[i].name << ";" << inps[i] << endl;
    }
    out.close();
    if (!before) {
        out.open((filename + "_out.csv").c_str(), ios::trunc);
        for (int i = 0; i < surrogateIO.output_names.size(); ++i) {
            double val = 101;
            if (i < x.size())
                val = x[i];
            out << remove_price_suffix(surrogateIO.output_names[i]) << ";" << val << endl;
        }
        out.close();
    }
} 

static void mkIDs() {
    int msz = marketdata.products.size();
    int isz = investdata.invests.size();
    for (int i = 0 ; i < msz; i++)
        marketId[marketdata.products[i].name]=i;
    for (int i = 0; i< isz; i++)
        investId[investdata.invests[i].name]=i;
    return;
}

double RegSurrogate::getValOfIndex(int id) {
	return x[id];
}

void RegSurrogate::setupSurrogate(RegGlobalsInfo* g){
    this->g = g;
    inputsMinMax = surrogateIO.inMinMax;
    inputs.resize(surrogateIO.inputlinks.size(),0);
    int s = reftypes.size();
    for (int i = 0; i < s; ++i) {
        refnumber[reftypes[i]] = i;
    }
    mkIDs();
    int sz = surrogateIO.output_names.size();
    for (int i = 0; i < sz; ++i) {
        pair<objtype, int> res = result_type(surrogateIO.output_names[i]);
        switch (res.first) {
        case UNKNOWN: g->Surrogate_Extra_outnames.push_back(surrogateIO.output_names[i]);
            break;
        default:;
        }
    }
    int dk = 3; //surrogate
    for (unsigned int i = 0;i < surrogateIO.inputlinks.size();i++) {
        RegLinkObject* link = readLink(i, dk);
        //input_links.push_back(link);
    }
 }

bool
RegSurrogate::updateInputValues() {
    bool changed=false;
    for (unsigned i=0;i<invest_links.size();i++)
        if (invest_links[i]->trigger()) changed=true;
    for (unsigned i=0;i<market_links.size();i++)
        if (market_links[i]->trigger()) changed=true;
    for (unsigned i=0;i<reference_links.size();i++)
        if (reference_links[i]->trigger()) changed=true;
    for (unsigned i=0;i<number_links.size();i++)
        if (number_links[i]->trigger()) changed=true;
	if (g->HAS_SOILSERVICE){
		for (unsigned i=0;i<yield_links.size();i++)
			if (yield_links[i]->trigger()) changed=true;
	}
	return changed;
}
void
RegSurrogate::updateLand() {
    // land_links is pointer to object of type RegLinkObject
    for (unsigned i=0;i < land_links.size();i++)
        land_links[i]->trigger();
    
    updateReferences();
}

void
RegSurrogate::updateReferences() {
    for (unsigned i = 0; i < reference_links.size(); i++)
        reference_links[i]->trigger();
}

void
RegSurrogate::updateYield() {
    for (unsigned i=0;i < yield_links.size();i++)
        yield_links[i]->trigger();
}

void RegSurrogate::updatePaymentEntitlement() {
    for (unsigned i=0;i < incomepay_links.size();i++)
        incomepay_links[i]->trigger();
}

void RegSurrogate::adjustInputs(vector<float>& inps) {
    //with inputsMinMax
    auto sz = inps.size();
    string name;
    float mi, ma;
    float v;
    for (auto i = 0; i < sz; ++i) {
        name = get<0>(inputsMinMax[i]);
        mi = get<1>(inputsMinMax[i]);
        ma = get<2>(inputsMinMax[i]);
        v = inps[i];
        if (v > ma ) {
            //cout << "WARNING "<<"("<< v << ">" << ma<< "): #" << name << " is larger than training maximum !, corrected to max \n";
            inps[i] = ma;
        }else if (v<mi) {
            //cout << "WARNING " << "(" << v << "<" << mi << "): #" << name << " is smaller than training minimum !, corrected to min \n";
            inps[i] = mi;
        }
    }

    /*
    if (inps[3]>2)  //chiselPlough
        inps[3]=2;
    if (inps[79] < 15.02)   //nArabLand
        inps[79] = 15.02;
    */
}

void RegSurrogate::adjustOutputs(vector<float>& x, const vector<float>& inp ) {
    x[40] = x[48] = x[49] = 0;  //solidManDist, potaStore500t, bunkerSilo450
    if (x[61] > 0) {   // EC_INTEREST
        x[61] /= inp[128];
    }
}

double RegSurrogate::LpSurrogate(RegProductList* PList, vector<int >& ninv, vector<double>& extras, int maxofffarmlu) {
    extras.clear();
    if (PList->setUsePriceExpectation(true)) { // hier war der Fehler
        for (unsigned i = 0;i < market_links.size();i++)
            market_links[i]->trigger();
    }

#ifndef NDEBUG
    bool abool0 = false;

    if (abool0) {
        stringstream ss;
        ss << g->tIter << "_" << g->tFarmName << "_" << g->tFarmId << "_debug0.txt";
        debug(ss.str().c_str(),true);
    }
#endif

    vector<float> inp;
    /*
    for (int i = 0;i < g->SurrogateParas.dim_in; ++i) {
        inp.push_back(1000);
    }
   //*/
 
    inp.resize(g->SurrogateParas.dim_in);
    int c = 0;
    for (auto v : inputs) {
        //inp.push_back(v);
        if (c < inp.size())
            inp[c++] = v;
    }
    adjustInputs(inp);

    //*/ 
    vector<float> outp;
    //outp.resize(g->SurrogateParas.dim_out);

    x.resize(g->SurrogateParas.dim_out);
    
    predict(inp, x, g->SurrogateParas.input_name, g->SurrogateParas.output_name);
    //if (g->tPhase == SimPhase::INVEST)
       //debugCSV(inp, "Farm_" + to_string(g->tFarmId), false);
    
    if (x[0] <= 0) {
        debugPredict("Farm_"+to_string(g->tFarmId),inp, x);
        cout << "GDB<=0" << endl; // exit(5);
    }

    adjustOutputs(x, inp);
   /*
    outp = x;
    std::cout << outp.size() << std::endl; 

    std::cout << setprecision(8);
    std::cout << std::defaultfloat << std::endl;
    std::cout << std::setw(30);
    int cnt = 0;
    for (int i = 0;i < g->SurrogateParas.dim_out;++i) {
        std::cout << std::setw(20) << outp[i];
        //++cnt;
        //if (cnt % 6 == 0) std::cout << std::endl; 
        //else std::cout << "\t";
    }
    //*/

#ifndef NDEBUG
	bool abool = false;
	if (abool) {
		stringstream ss;
		ss << g->tIter << "_" << g->tFarmName << "_" << g->tFarmId << "_debug.txt";
		debug(ss.str().c_str());
	}
#endif

        objtype otype;
        double td; long long int ti;
        int sz = surrogateIO.output_names.size();
     for (int i=0; i<sz; ++i){
        pair<objtype,int> res = result_type(surrogateIO.output_names[i]);
        switch (res.first) {
        case PROD: 
            //  PList->setUnitsProducedOfNumber(s,x[d]);
			td = x[i];
			ti = (long long int)(td * 1E+7 + 0.5);
			td = ti * 1E-7;
			PList->setUnitsProducedOfNumber(res.second, td);
			break;
		case INVEST: // number of investments of investment number
			ninv[res.second] = (int)(x[i] + 0.5);
			break;
        case OBJVALUE:
            objval = x[i];
            break;
        case UNKNOWN:
            extras.push_back(x[i]);
            break;
		default:
			// for debug 
            //cout << "ELPSURROGATE: solution assignment!" << endl;
            //exit(3);
            ;
		}
	}
    
#ifndef NDEBUG1
	if (g->DebMip)	{
		bool mipcond = false;
		bool itercond = g->uIter==-1 || g->uIter == g->tIter;
		bool idcond = g->uFarmId==-1 || g->uFarmId == g->tFarmId;
		bool namecond = g->uFarmName.compare("ALL")==0 || g->uFarmName.compare(g->tFarmName)==0;
		bool phasecond =g->tPhase!=SimPhase::BETWEEN && ( g->uPhase == SimPhase::ALL || g->uPhase == g->tPhase);

		mipcond = itercond && idcond && namecond && phasecond;
		string phasestr;
		switch (g->tPhase) {
		case SimPhase::FUTURE: phasestr = "future"; break;
		case SimPhase::INVEST: phasestr = "invest"; break;
		case SimPhase::LAND: phasestr = "land"; break;
		case SimPhase::PRODUCT: phasestr = "product"; break;
		default:phasestr = "impossible";
		}

		if ( mipcond ) {
			string indstr;
			if (g->tPhase == SimPhase::LAND ){
				int r = g->mapMIP[tuple(g->tIter, g->tFarmId, SimPhase::LAND)]++;
				indstr = to_string(r);// g->tInd_land);
				//++g->tInd_land;
			}else if (g->tPhase == SimPhase::FUTURE) {
				int r = g->mapMIP[tuple(g->tIter, g->tFarmId, SimPhase::FUTURE)]++;
				indstr = to_string(r);// g->tInd_future);
				//++g->tInd_future;
			}

			stringstream ts;
			string debdir = "DebMIPs\\";
		
			string tinputdir;
			tinputdir=rtrim(inputdir, '\\');
				
			fs::path inpdir(tinputdir);
			fs::path pdir = inpdir.parent_path();
			string dstr = pdir.string() + "\\" + debdir;
			fs::path ddir(dstr);
			if (!fs::exists(ddir))
				fs::create_directories(ddir);
			ts << dstr;
			
			ts << "It_" << g->tIter << "_Id_" << g->tFarmId << "_"
				<< g->tFarmName << "_" << phasestr;
			if (indstr.length()>0)
				ts << "_" << indstr << ".txt";
			else
				ts << ".txt";

			debug(ts.str().c_str());
		}
	}
#endif

	if (g->SDEBUG1 || g->SDEBUG2) {
		stringstream sstr;
		sstr << g->tInd;
		string str = sstr.str();
		if (g->SDEBUG1)
			debug("debug-r" + str + ".txt");
		if (g->SDEBUG2)
			debug("debug-n" + str + ".txt");
	}

	bool deb = false;
	
	if (deb)
        debug("ldebug2.txt");
    
    return objval;
} 

RegLinkObject* RegSurrogate::mklink(onelink& lk, int dn, int dk) {
    int vk;
    string skString= lk.linktype;
    std::transform(skString.begin(), skString.end(), skString.begin(),
               (int(*)(int)) std::toupper);
    if (skString=="MARKET") {
		int sn = marketId[lk.numbertype];
        if (lk.numbertype == "hourly_wage__mean")
            sn = marketId["V_HIRED_LABOUR"];
        string vkString= lk.valuetype;
        std::transform(vkString.begin(), vkString.end(), vkString.begin(),
               (int(*)(int)) std::toupper);
        
        if (vkString=="C")
            vk=0;
        if (vkString=="P")
            vk=1;
        if (vkString=="PE")
            vk=2;
        if (vkString=="GM")
            vk=3;
    //TODO
        if (vkString == "NNP")
            vk = 4;
        if (vkString == "NNY")
            vk = 5; 

		double f = lk.factor;
        RegLinkMarketObject *link = new RegLinkMarketObject(dn,dk,sn,vk,f);
        market_links.push_back(link);
        return link;
    }
    if (skString=="INVEST") {
        int sn;
        if (dk == 1) //capacity links
           sn = atoi(lk.numbertype.c_str());
        else {
			sn = investId[lk.numbertype];
		}
		string vkString= lk.valuetype;
        std::transform(vkString.begin(), vkString.end(), vkString.begin(),
               (int(*)(int)) std::toupper);
        if (vkString=="CAP")
            vk=0;
        if (vkString=="LE")
            vk=1;
        if (vkString=="BE")
            vk=2;
        if (vkString=="AC")
            vk=3;
        if (vkString=="NORMCAP")
            vk=4;
        //TODO
        if (vkString == "NNAGE")
            vk = 5;  
        if (vkString == "NNCAP")
            vk = 6;
        double f = lk.factor;
        RegLinkInvestObject *link =  new RegLinkInvestObject(dn,dk,sn,vk,f);
        invest_links.push_back(link);
        return link;
    }
    
    if (skString=="REFERENCE") {
        string vkString= lk.valuetype;
        std::transform(vkString.begin(), vkString.end(), vkString.begin(),
               (int(*)(int)) std::toupper);
        vk = refnumber[vkString];
        
		double f = lk.factor;
        RegLinkReferenceObject *link =   new RegLinkReferenceObject(dn,dk,vk,f);
        reference_links.push_back(link);
        if (vkString.compare("INCOMEPAY")==0 || vkString.compare("UNMODINCOMEPAY")==0)
            incomepay_links.push_back(link);
        return link;
    }
    if (skString=="NUMBER") {
        double v = atof(lk.numbertype.c_str());
        RegLinkNumberObject *link =  new RegLinkNumberObject(dn,dk,v);
        number_links.push_back(link);
        return link;
    }
    if (skString=="LAND") {
        string vkString= lk.valuetype;
        int i=0;
        for (i=0;i<g->NO_OF_SOIL_TYPES;i++)  {
            if (vkString==g->NAMES_OF_SOIL_TYPES[i])
                break;
        }
        double f = lk.factor;
        RegLinkLandObject *link = new RegLinkLandObject(dn,dk,i,f);
        land_links.push_back(link);
        return link;
    }

  if (g->HAS_SOILSERVICE) {
	if (skString=="YIELD") {
		string vkString= lk.valuetype;
		std::transform(vkString.begin(), vkString.end(), vkString.begin(),
               (int(*)(int)) std::toupper);
        int sn = marketId[lk.numbertype];

		if (vkString=="SY")
            vk=0;
        if (vkString=="P")
            vk=1;
        if (vkString=="K")
            vk=2;
        if (vkString=="PE")
            vk=3;
		if (vkString=="ENV")
            vk=4;
		if (vkString=="SN")
			vk=5;
		
        double f = lk.factor;
        RegLinkYieldObject *link = new RegLinkYieldObject(dn,dk,sn,vk,f);
        yield_links.push_back(link);
        return link;
    }
  }
    return NULL;
}

RegLinkObject* RegSurrogate::readLink(int dn,int dk) {
    onelink lk;
    switch (dk) {
    case 3:     // surrogate links
                lk = surrogateIO.inputlinks[dn];
                break;
    default: 
        break;
    }

    return mklink(lk,dn,dk);
}

RegSurrogate* RegSurrogate::create() {
    return new RegSurrogate();
}

RegSurrogate* RegSurrogate::clone() {
    RegSurrogate* cl=new RegSurrogate(*this);
    cl->flat_copy=true;
    return cl;
}

RegSurrogate* RegSurrogate::clone(RegGlobalsInfo* G) {
    RegSurrogate* n=create();
    (*n).g=G;
    (*n).refnumber = refnumber;
    (*n).inputs.resize(inputs.size());
    (*n).inputsMinMax = inputsMinMax;

    for (unsigned i=0;i<invest_links.size();i++) {
        RegLinkInvestObject *link=new RegLinkInvestObject(*invest_links[i]);
        (*n).invest_links.push_back(link);
        input_links.push_back(link);
    }
    for (unsigned i=0;i<market_links.size();i++) {
        RegLinkMarketObject *link=new RegLinkMarketObject(*market_links[i]);
        (*n).market_links.push_back(link);
        input_links.push_back(link);
    }

	//soil service yield links 
	for (unsigned i=0;i<yield_links.size();i++) {
        RegLinkYieldObject *link=new RegLinkYieldObject(*yield_links[i]);
        (*n).yield_links.push_back(link);
        input_links.push_back(link);
    }
    for (unsigned i=0;i<land_links.size();i++) {
        RegLinkLandObject *link=new RegLinkLandObject(*land_links[i]);
        (*n).land_links.push_back(link);
        input_links.push_back(link);
    }
    for (unsigned i=0;i<reference_links.size();i++) {
        RegLinkReferenceObject *link=new RegLinkReferenceObject(*reference_links[i]);
        (*n).reference_links.push_back(link);
        input_links.push_back(link);
        if (link->getSourceNumber()==7 || link->getSourceNumber()==8)
            (*n).incomepay_links.push_back(link);
    }
    for (unsigned i=0;i<number_links.size();i++) {
        RegLinkNumberObject *link=new RegLinkNumberObject(*number_links[i]);
        (*n).number_links.push_back(link);
        input_links.push_back(link);
    }
    return n;
}

void
RegSurrogate::initLinks(RegFarmInfo*f,
                     RegInvestList* il,
                     RegProductList* pl,
                     RegLabourInfo *lab ) {
	farm = f;
    
    for (unsigned i=0;i<invest_links.size();i++) {
        invest_links[i]->init(il, &(*inputs.begin()));
    }
    
    for (unsigned i=0;i<market_links.size();i++) {
        market_links[i]->init(pl, &(*inputs.begin()));
    }
    for (unsigned i=0;i<land_links.size();i++) {
        land_links[i]->init(f, &(*inputs.begin()));
    }

	//soil service
	for (unsigned i=0;i<yield_links.size();i++) {
        yield_links[i]->init(f, &(*inputs.begin()));
    }

    //  "Management_Coeff", "LIQUIDITY", "LAND", "LABOUR", "FINANCIALRULE", "INCOMEPAY",
    //  "UNMODINCOMEPAY", "TRANCH1WIDTH", "TRANCH2WIDTH", "TRANCH3WIDTH", "TRANCH4WIDTH",
    //  "TRANCH5WIDTH", "TRANCH1DEG", "TRANCH2DEG", "TRANCH3DEG", "TRANCH4DEG", 
    //  "TRANCH5DEG", "LABOUR_SUBSTITUTION", "FAMILYLABOUR", "HIREDLABOUR",  "EC_INTEREST",
    // "ARABLENORMAL", "ARABLEREDZONE"
    
    refsources = { &(f->management_coefficient), &(f->liquidity), &(f->land_input), &lab->labour_capacity, &(f->financing_rule), &(f->modulated_income_payment),
                  &(f->income_payment_farm), &(g->TRANCH_1_WIDTH), &(g->TRANCH_2_WIDTH), &(g->TRANCH_3_WIDTH), &(g->TRANCH_4_WIDTH),
                  &(g->TRANCH_5_WIDTH), &(g->TRANCH_1_DEG), &(g->TRANCH_2_DEG), &(g->TRANCH_3_DEG), &(g->TRANCH_4_DEG),
                  &(g->TRANCH_5_DEG), &il->labSubstitution, &lab->family_labour, &lab->fix_onfarm_labour, 
                  &(g->EQ_INTEREST), &(f->arable_redzone), &(f->arable_nonRedzone)};

    for (unsigned i=0;i<reference_links.size();i++) {
        double* source;
        int num = reference_links[i]->getSourceNumber();
        source = refsources[num];
        	
        reference_links[i]->init(source,&(*inputs.begin()));
    }
    for (unsigned i=0;i<number_links.size();i++) {
        number_links[i]->init(&(*inputs.begin()));
    }
}
