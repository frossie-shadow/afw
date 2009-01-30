// -*- lsst-c++ -*-
//
//##====----------------                                ----------------====##/
//
//! \file
//! \brief  Testing of IO via the persistence framework for Source and SourceVector.
//
//##====----------------                                ----------------====##/

#include <sys/time.h>
#include <iostream>
#include <sstream>
#include <cstring>
#include <stdexcept>


#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE SourceIO

#include "boost/test/unit_test.hpp"

#include "lsst/daf/base.h"
#include "lsst/daf/persistence.h"
#include "lsst/pex/exceptions.h"
#include "lsst/pex/policy/Policy.h"

#include "lsst/afw/detection/Source.h"
#include "lsst/afw/formatters/Utils.h"

using lsst::daf::base::PropertySet;
using lsst::daf::base::Persistable;
using lsst::daf::persistence::LogicalLocation;
using lsst::daf::persistence::Persistence;
using lsst::daf::persistence::Storage;
using lsst::pex::policy::Policy;

namespace afwFormatters = lsst::afw::formatters;
using namespace lsst::afw::detection;


#define Assert(pred, msg) do { if (!(pred)) { throw LSST_EXCEPT(lsst::pex::exceptions::RuntimeErrorException, (msg)); } } while (false)


static std::string const makeTempFile() {
    char name[32];
    std::strncpy(name, "Source_XXXXXX", 31);
    name[31] = 0;
    int const fd = ::mkstemp(name);
    Assert(fd != -1, "Failed to create temporary file");
    ::close(fd);
    return std::string(name);
}


static void initTestData(SourceVector & v, int sliceId = 0) {
    v.clear();
    v.reserve(NUM_SOURCE_NULLABLE_FIELDS + 2);
    Source data;
    for (int i = 0; i != NUM_SOURCE_NULLABLE_FIELDS + 2; ++i) {
        // make sure each field has a different value, and that IO for each nullable field is tested
        // Note: Source ids are generated in ascending order
        int j = i*NUM_SOURCE_NULLABLE_FIELDS;
        data.setSourceId(j + sliceId*(NUM_SOURCE_NULLABLE_FIELDS + 2)*64 + 1);
        data.setAmpExposureId(static_cast<int64_t>(j +  2));
        data.setFilterId(-1);
        data.setObjectId(static_cast<int64_t>(j +  4));
        data.setMovingObjectId(static_cast<int64_t>(j +  5));
        data.setProcHistoryId(-1);
        data.setRa(static_cast<double>(j +  8));
        data.setRaErrForDetection(static_cast<float>(j +  9));
        data.setRaErrForWcs(static_cast<float>(j + 10));
        data.setDec(static_cast<double>(j + 11));
        data.setDecErrForDetection(static_cast<float>(j + 12));
        data.setDecErrForWcs(static_cast<float>(j + 13));
        data.setXFlux(static_cast<double>(j + 14));
        data.setXFluxErr(static_cast<double>(j + 15));
        data.setYFlux(static_cast<double>(j + 16));
        data.setYFluxErr(static_cast<double>(j + 17));
        data.setRaFlux(static_cast<double>(j + 18));
        data.setRaFluxErr(static_cast<double>(j + 19));
        data.setDecFlux(static_cast<double>(j + 20));
        data.setDecFluxErr(static_cast<double>(j + 21));
        data.setXPeak(static_cast<double>(j + 22));
        data.setYPeak(static_cast<double>(j + 23));
        data.setRaPeak(static_cast<double>(j + 24));
        data.setDecPeak(static_cast<double>(j + 25));
        data.setXAstrom(static_cast<double>(j + 26));
        data.setXAstromErr(static_cast<double>(j + 27));
        data.setYAstrom(static_cast<double>(j + 28));
        data.setYAstromErr(static_cast<double>(j + 29));        
        data.setRaAstrom(static_cast<double>(j + 30));
        data.setRaAstromErr(static_cast<double>(j + 31));
        data.setDecAstrom(static_cast<double>(j + 32));
        data.setDecAstromErr(static_cast<double>(j + 33));                
        data.setTaiMidPoint(static_cast<double>(j + 34));
        data.setTaiRange(static_cast<float>(j + 35));
        data.setFwhmA(static_cast<float>(j + 36));
        data.setFwhmB(static_cast<float>(j + 37));
        data.setFwhmTheta(static_cast<float>(j + 38));       
        data.setPsfMag(static_cast<double>(j + 39));
        data.setPsfMagErr(static_cast<float>(j + 40));
        data.setApMag(static_cast<double>(j + 41));
        data.setApMagErr(static_cast<float>(j + 42));
        data.setModelMag(static_cast<double>(j + 43));
        data.setModelMagErr(static_cast<float>(j + 44));
        data.setPetroMag(static_cast<double>(j + 45));
        data.setPetroMagErr(static_cast<float>(j + 46));
        data.setInstMag(static_cast<double>(j + 47));
        data.setInstMagErr(static_cast<double>(j + 48));
        data.setNonGrayCorrMag(static_cast<double>(j + 49));
        data.setNonGrayCorrMagErr(static_cast<double>(j + 50));
        data.setAtmCorrMag(static_cast<double>(j + 51));
        data.setAtmCorrMagErr(static_cast<double>(j + 52));
        data.setApDia(static_cast<float>(j + 53));
        data.setSnr(static_cast<float>(j + 54));
        data.setChi2(static_cast<float>(j + 55));
        data.setSky(static_cast<float>(j + 56));
        data.setSkyErr(static_cast<float>(j + 57));
        data.setFlagForAssociation(1);
        data.setFlagForDetection(2);
        data.setFlagForWcs(3);

        if (i < NUM_SOURCE_NULLABLE_FIELDS) {
            data.setNotNull();
            data.setNull(i);
        } else if ((i & 1) == 0) {
            data.setNotNull();
        } else {
            data.setNull();
        }
        
        Source::Ptr sourcePtr(new Source(data));
        v.push_back(sourcePtr);
    }
}

// Make at least a token attempt at generating a unique visit id
// (in-db table name collisions could cause spurious testcase failures)
static int createVisitId() {
    struct timeval tv;
    ::gettimeofday(&tv, 0);
    return abs(static_cast<int>(tv.tv_sec));
}

static PropertySet::Ptr createDbTestProps(
    int const sliceId,
    int const numSlices,
    std::string const & itemName
) {
    Assert(sliceId < numSlices && numSlices > 0, "invalid slice parameters");

    PropertySet::Ptr props(new PropertySet);

    if (numSlices > 1) {
        props->add("Source.isPerSliceTable", true);
        props->add("Source.numSlices", numSlices);
    }
    int visitId = createVisitId();
    props->add("visitId",  visitId);
    props->add("exposureId", static_cast<int64_t>(visitId*2));
    props->add("ccdId", static_cast<int64_t>(5));
    props->add("universeSize", numSlices);
    props->add("sliceId",  sliceId);
    props->add("itemName", itemName);
    return props;
}

static void testBoost(void) {
    // Create a blank Policy and PropertySet.
    Policy::Ptr      policy(new Policy);
    PropertySet::Ptr props = createDbTestProps(0,1,"Source");

    // Setup test location
    LogicalLocation loc(makeTempFile());

    // Intialize test data
    SourceVector dsv;
    initTestData(dsv);
    PersistableSourceVector::Ptr persistPtr(new PersistableSourceVector(dsv));    
    Persistence::Ptr pers = Persistence::getPersistence(policy);

    // write out data
    {
        Storage::List storageList;
        storageList.push_back(pers->getPersistStorage("BoostStorage", loc));
        pers->persist(*persistPtr, storageList, props);
    }

    // read in data
    {
        Storage::List storageList;
        storageList.push_back(pers->getRetrieveStorage("BoostStorage", loc));
        Persistable::Ptr p = pers->retrieve("PersistableSourceVector", storageList, props);
        BOOST_CHECK_MESSAGE(p.get() != 0, "Failed to retrieve Persistable");
        PersistableSourceVector::Ptr persistVec =
            boost::dynamic_pointer_cast<PersistableSourceVector, Persistable>(p);
        BOOST_CHECK_MESSAGE(persistVec.get() != 0, "Couldn't cast to PersistableSourceVector");
        BOOST_CHECK_MESSAGE(*persistVec == dsv, 
            "persist()/retrieve() resulted in PersistableSourceVector corruption");
    }
    ::unlink(loc.locString().c_str());
}





// comparison operator used to sort Source in id order
struct SourceLessThan {
    bool operator()(Source::Ptr const & d1, Source::Ptr const & d2) {
        return d1->getId() < d2->getId();
    }
};


static void testDb(std::string const & storageType) {
    // Create the required Policy and PropertySet
    Policy::Ptr policy(new Policy);
    PropertySet::Ptr props = createDbTestProps(0, 1, "Source");
    
    Persistence::Ptr pers = Persistence::getPersistence(policy);
    LogicalLocation loc("mysql://lsst10.ncsa.uiuc.edu:3306/source_test");
    
    // 1. Test on a single Source
    Source::Ptr ds(new Source());
    ds->setId(2);
    SourceVector dsv;
    dsv.push_back(ds);
    PersistableSourceVector::Ptr persistPtr(new PersistableSourceVector(dsv));
    // write out data
    {
        Storage::List storageList;
        storageList.push_back(pers->getPersistStorage(storageType, loc));
        pers->persist(*persistPtr, storageList, props);
    }
    // and read it back in (in a SourceVector)
    {
        Storage::List storageList;
        storageList.push_back(pers->getRetrieveStorage(storageType, loc));
        Persistable::Ptr p = pers->retrieve("PersistableSourceVector", storageList, props);
        BOOST_CHECK_MESSAGE(p.get() != 0, "Failed to retrieve Persistable");
        PersistableSourceVector::Ptr persistVec = 
            boost::dynamic_pointer_cast<PersistableSourceVector, Persistable>(p);
        BOOST_CHECK_MESSAGE(persistVec.get() != 0, "Couldn't cast to PersistableSourceVector");
        SourceVector vec = persistVec->getSources();
        BOOST_CHECK_MESSAGE(*vec.at(0) == *dsv[0], 
            "persist()/retrieve() resulted in PersistableSourceVector corruption");
    }
    afwFormatters::dropAllVisitSliceTables(loc, policy, props);

    // 2. Test on a SourceVector
    dsv.clear();
    initTestData(dsv);
    persistPtr->setSources(dsv);
    // write out data
    {
        Storage::List storageList;
        storageList.push_back(pers->getPersistStorage(storageType, loc));
        pers->persist(*persistPtr, storageList, props);
    }
    // and read it back in
    {
        Storage::List storageList;
        storageList.push_back(pers->getRetrieveStorage(storageType, loc));
        Persistable::Ptr p = pers->retrieve("PersistableSourceVector", storageList, props);
        BOOST_CHECK_MESSAGE(p.get() != 0, "Failed to retrieve Persistable");
        PersistableSourceVector::Ptr persistVec = 
                boost::dynamic_pointer_cast<PersistableSourceVector, Persistable>(p);
        BOOST_CHECK_MESSAGE(persistVec.get() != 0, "Couldn't cast to PersistableSourceVector");
        SourceVector vec(persistVec->getSources());
        
        // sort in ascending id order (database does not give any ordering guarantees
        // in the absence of an ORDER BY clause)
        std::sort(vec.begin(), vec.end(), SourceLessThan());
        BOOST_CHECK_MESSAGE(vec.size() == dsv.size(),
            "persist()/retrieve() resulted in PersistableSourceVector corruption");
    
        for (size_t i =0; i<vec.size();i++){
            if (*vec[i] != *dsv[i]){
                BOOST_ERROR("persist()/retrieve() resulted in PersistableSourceVector corruption");
                break;
            }
        }
    }
    afwFormatters::dropAllVisitSliceTables(loc, policy, props);
}


static void testDb2(std::string const & storageType) {
    // Create the required Policy and PropertySet
    Policy::Ptr policy(new Policy);
    std::string policyRoot("Formatter.PersistableSourceVector");
    // use custom table name patterns for this test
    policy->set(policyRoot + ".Source.perVisitTableNamePattern", "tmp_test_Source_%1%");
    policy->set(policyRoot + ".Source.perSliceAndVisitTableNamePattern", "tmp_test_Source_%1%_%2%");

    Policy::Ptr nested(policy->getPolicy(policyRoot));

    Persistence::Ptr pers = Persistence::getPersistence(policy);
    LogicalLocation loc("mysql://lsst10.ncsa.uiuc.edu:3306/source_test");

    SourceVector all;
    int const numSlices = 3;
    PropertySet::Ptr props = createDbTestProps(0, numSlices, "Source");

    // 1. Write out each slice table seperately
    for (int sliceId = 0; sliceId < numSlices; ++sliceId) {
        props->set("sliceId", sliceId);
        SourceVector dsv;
        initTestData(dsv, sliceId);

        all.insert(all.end(), dsv.begin(), dsv.end());
        PersistableSourceVector persistVec;
        persistVec.setSources(dsv);
        Storage::List storageList;
        storageList.push_back(pers->getPersistStorage(storageType, loc));
        pers->persist(persistVec, storageList, props);
    }

    // 2. Read in all slice tables - simulates association pipeline
    //    gathering the results of numSlices image processing pipeline slices
    Storage::List storageList;
    storageList.push_back(pers->getRetrieveStorage(storageType, loc));
    Persistable::Ptr p = pers->retrieve("PersistableSourceVector", storageList, props);
    BOOST_CHECK_MESSAGE(p.get() != 0, "Failed to retrieve Persistable");
    PersistableSourceVector::Ptr persistPtr = 
        boost::dynamic_pointer_cast<PersistableSourceVector, Persistable>(p);
    BOOST_CHECK_MESSAGE(persistPtr.get() != 0, "Couldn't cast to PersistableSourceVector");
    
    // sort in ascending id order (database does not give any ordering guarantees
    // in the absence of an ORDER BY clause)
    SourceVector vec = persistPtr->getSources();
    std::sort(vec.begin(), vec.end(), SourceLessThan());
    BOOST_CHECK_MESSAGE(vec.size() == all.size(),
        "persist()/retrieve() resulted in PersistableSourceVector corruption");

    for (size_t i =0; i< vec.size();i++){
        if (*vec[i] != *all[i]){
            BOOST_ERROR("persist()/retrieve() resulted in PersistableSourceVector corruption");
            break;
        }            
    }    
    afwFormatters::dropAllVisitSliceTables(loc, nested, props);
}

BOOST_AUTO_TEST_CASE(SourceEquality) {
    Source::Ptr a(new Source), b(new Source);
    a->setId(3);
    BOOST_CHECK_MESSAGE(*a != *b && *b != *a, "field equality fails");
    b->setId(3);
    BOOST_CHECK_MESSAGE(*a == *b && *b == *a, "field equality fails");
    a->setNotNull(lsst::afw::detection::MOVING_OBJECT_ID);
    BOOST_CHECK_MESSAGE(*a != *b && *b != *a, "field equality fails");    
    a->setMovingObjectId(5);
    BOOST_CHECK_MESSAGE(*a != *b && *b != *a, "field equality fails");    
    b->setNotNull(lsst::afw::detection::MOVING_OBJECT_ID);
    BOOST_CHECK_MESSAGE(*a != *b && *b != *a, "field equality fails");    
    b->setMovingObjectId(5);    
    BOOST_CHECK_MESSAGE(*a == *b && *b == *a, "field equality fails");
    
    SourceVector av, bv;
    av.push_back(a);
    PersistableSourceVector apv(av);
    BOOST_CHECK(apv.getSources()[0]->isNull(0) == a->isNull(0));
    BOOST_CHECK(apv.getSources()[0]->isNull(1) == a->isNull(1));
    bv.push_back(b);
    
}

BOOST_AUTO_TEST_CASE(SourceIO) {
    try {
        testBoost();
        if (lsst::daf::persistence::DbAuth::available()) {
            testDb("DbStorage");
            testDb("DbTsvStorage");
            testDb2("DbStorage");
            testDb2("DbTsvStorage");
        }
        BOOST_CHECK_MESSAGE(lsst::daf::base::Citizen::census(0) == 0, "Detected memory leaks");
    } catch(std::exception const & ex) {
        BOOST_FAIL(ex.what());
    }
}
