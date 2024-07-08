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

//---------------------------------------------------------------------------
// RegSurrogate.h
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
#ifndef RegSurrogateH
#define RegSurrogateH
//---------------------------------------------------------------------------
#include <vector>
#include <stdlib.h>
#include <fstream>
#include "RegGlobals.h"
#include "RegProduct.h"
#include "RegStructure.h"
#include "RegInvest.h"
#include "RegLink.h"

#include "textinput.h"

#include "RegLabour.h"

//---------------------------------------------------------------------------

/** RegSurrogate class.
    @short class that defines the functions related to the surrogate model
 */

class RegLinkInvestObject;
class RegLinkMarketObject;

class RegRegionInfo;
class RegProductList;

class RegLinkYieldObject;

class RegLinkReferenceObject;
class RegLinkNumberObject;
class RegLinkLandObject;
class RegLinkObject;
class RegFarmInfo;

class RegSurrogate {
protected:
	RegFarmInfo* farm;
	
    //inputs for surrogate model
    vector<double> inputs;

    bool flat_copy;
    /// pointer to globals
    RegGlobalsInfo* g;
       
    // status of the solution
    long     stat;
    // objective function value
    double   objval;
        
    // array that contains the optimal values of the primal variables
    vector<float> x;
   
    // INITIAL LINK LISTS
    // they are necessary for initialisation, but not thereafter
    vector<RegLinkInvestObject*> invest_links;
    vector<RegLinkMarketObject*> market_links;
    vector<RegLinkReferenceObject*> reference_links;
    const vector<string> reftypes = { "MANAGEMENT", "LIQUIDITY", "LAND", "LABOUR", "FINANCIALRULE", "INCOMEPAY",
        "UNMODINCOMEPAY", "TRANCH1WIDTH", "TRANCH2WIDTH",  "TRANCH3WIDTH", "TRANCH4WIDTH",
        "TRANCH5WIDTH", "TRANCH1DEG", "TRANCH2DEG", "TRANCH3DEG", "TRANCH4DEG",
        "TRANCH5DEG", "LABOUR_SUBSTITUTION", "FAMILYLABOUR", "HIREDLABOUR", "EC_INTEREST",
        "ARABLENORMAL", "ARABLEREDZONE"};
    map <string, int> refnumber;
    vector<double*> refsources;
    vector<RegLinkNumberObject*> number_links;
    vector<RegLinkLandObject *> land_links;
    
    //soil service yield links
    vector<RegLinkYieldObject*> yield_links;
    vector<RegLinkObject *> incomepay_links;

    // Lists containing the destinations of values
    vector<RegLinkObject*> input_links;

      
    /** method that connects sources to destinations;
        values are retrieved from the input file and written 
        at the appropriate position
        @param dn destination
        @param dk destination kind
   */
    RegLinkObject* readLink(int dn, int dk);
    
    //gemeinsames 
    RegLinkObject* mklink(onelink& lk, int dn, int dk);

    RegSurrogate* obj_backup;
public:
    void    setupSurrogate(RegGlobalsInfo* G);
    double getValOfIndex(int);

    /** Surrogate model predict method
        @param PList pointer to FarmProductList in RegFarmInfo
        @param reference to inum_vector in RegFarmInfo, which keeps track
                of the number of investments in to one object
        @return objective function value
    */
    virtual double  LpSurrogate(RegProductList* PList,vector<int >& ninv, vector<double>& extras, int famlabour);

    /** input variables change during runtime, actual values are retrieved   */
    bool    updateInputValues();

    void debug(string,bool);
    void debugCSV(vector<float>, string, bool);
    void initLinks(RegFarmInfo*,RegInvestList*,RegProductList*,RegLabourInfo*);
    /// to update only the land links
    void updateLand();
    void updateReferences();

	///update yield --soil service
	void updateYield();

    /// to update only the incomepay links
    void updatePaymentEntitlement();

    RegSurrogate();
    virtual RegSurrogate*  clone();

    virtual RegSurrogate*  clone(RegGlobalsInfo*);
    virtual RegSurrogate*  create();
    void setFlatCopy() {
        flat_copy=true;
    }
    /// Destructor
    ~RegSurrogate();
};
#endif
