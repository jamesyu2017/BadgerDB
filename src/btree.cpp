/**
 * @author See Contributors.txt for code contributors and overview of BadgerDB.
 *
 * @section LICENSE
 * Copyright (c) 2012 Database Group, Computer Sciences Department, University of Wisconsin-Madison.
 */

#include "btree.h"
#include "filescan.h"
#include "exceptions/bad_index_info_exception.h"
#include "exceptions/bad_opcodes_exception.h"
#include "exceptions/bad_scanrange_exception.h"
#include "exceptions/no_such_key_found_exception.h"
#include "exceptions/scan_not_initialized_exception.h"
#include "exceptions/index_scan_completed_exception.h"
#include "exceptions/file_not_found_exception.h"
#include "exceptions/end_of_file_exception.h"
#include "exceptions/file_exists_exception.h"
#include "exceptions/hash_not_found_exception.h"
#include "exceptions/page_not_pinned_exception.h"

//#define DEBUG

namespace badgerdb
{
// implement template to construct
template <class T>
// define copy
void copy(T& t1, T& t2){
    t1 = t2;
}

// specify char[] copy
template <>
void copy<char[STRINGSIZE]>(char (&t1)[STRINGSIZE], char (&t2)[STRINGSIZE]){
    strncpy(t1, t2, STRINGSIZE);
}

/*discuss comparison of equalto, smallerthan, biggerthan*/

// Type T equal comparison
template <class T>
bool equalTo(T t1, T t2){
    return t1 == t2;
}

// specialization char[STRINGSIZE] equal to comparison
template <>
bool equalTo<char[STRINGSIZE]>(char t1[STRINGSIZE], char t2[STRINGSIZE]){
    return strncmp(t1, t2, STRINGSIZE) == 0;
}

//  Type T bigger comparison
template <class T>
bool biggerThan(T t1, T t2){
    return t1 > t2;
}

// specialization char[STRINGSIZE] bigger comparison
template <>
bool biggerThan<char[STRINGSIZE]>(char t1[STRINGSIZE], char t2[STRINGSIZE]){
    return strncmp(t1, t2, STRINGSIZE) > 0;
}

// Type T small than comparison
template <class T>
bool smallerThan(T t1, T t2){
    return t1 < t2;
}

// specialiazation char[STRINGSIZE] smaller comparison
template <>
bool smallerThan<char[STRINGSIZE]>(char t1[STRINGSIZE], char t2[STRINGSIZE]){
    return strncmp(t1, t2, STRINGSIZE) < 0;
}

// -----------------------------------------------------------------------------
// BTreeIndex::BTreeIndex -- Constructor
// -----------------------------------------------------------------------------
/* The constructor checks 
 * if the specified index file exists.
 * and if index file name is constructed by concatenating the relational name with the offset of the attribute over which 
 * the index is built
 * */
BTreeIndex::BTreeIndex(const std::string & relationName,
		std::string & outIndexName,
		BufMgr *bufMgrIn,
		const int attrByteOffset,
		const Datatype attrType)
{
    // initialize all the local variables
    std::ostringstream idxStr;
    idxStr<<relationName<<'.'<<attrByteOffset;
    outIndexName = idxStr.str();
    headerPageNum = 1;
    bufMgr = bufMgrIn;
    attributeType = attrType;
    BTreeIndex::attrByteOffset = (int)attrByteOffset;
    
    if(attributeType == INTEGER){
        leafOccupancy = INTARRAYLEAFSIZE;
        nodeOccupancy = INTARRAYNONLEAFSIZE;
    }
    else if(attributeType == DOUBLE){
        leafOccupancy = DOUBLEARRAYLEAFSIZE;
        nodeOccupancy = DOUBLEARRAYNONLEAFSIZE;
    }
    else if(attributeType == STRING){
        leafOccupancy = STRINGARRAYLEAFSIZE;
        nodeOccupancy = STRINGARRAYNONLEAFSIZE;
    }
    else{
        std::cout << "Attribute Type does not fit the enum! \n";
        return;
    }
    scanExecuting = false;
    std::cout<< outIndexName << "\n";

    try{
        BlobFile* newFile = new BlobFile(outIndexName, true);
        file = (File*)newFile;
        Page* empty_page;
        PageId meta_pageId;
        bufMgr->allocPage(file, meta_pageId, empty_page);
        bufMgr->unPinPage(file, meta_pageId, true);
        bufMgr->readPage(file, meta_pageId, empty_page);
        IndexMetaInfo* metapage = (IndexMetaInfo*)empty_page;
        
        std::copy(relationName.begin(), relationName.end(), (*metapage).relationName);
        (*metapage).attrByteOffset = attrByteOffset;
        (*metapage).attrType = attrType;
        Page* root_page;
        bufMgr->allocPage(file, rootPageNum, root_page);
        bufMgr->unPinPage(file,rootPageNum, true);
        bufMgr->readPage(file, rootPageNum, root_page);
        (*metapage).rootPageNo = rootPageNum;
        (*metapage).isRootLeafPage = true;

        // create an empty leaf node at root
        if((*metapage).attrType == INTEGER){
            LeafNodeInt* leafNode = (LeafNodeInt*) root_page;
            (*leafNode).validEntryNo = 0;
        }
        else if((*metapage).attrType == DOUBLE){
            LeafNodeDouble* leafNode = (LeafNodeDouble*) root_page;
            (*leafNode).validEntryNo = 0;
        }
        else if((*metapage).attrType == STRING){
            LeafNodeString* leafNode = (LeafNodeString*) root_page;
            (*leafNode).validEntryNo = 0;
        }
        else{
            std::cout <<"INVALID type in constructor\n";
        }
        

        bufMgr->unPinPage(file, rootPageNum, true );
        bufMgr->unPinPage(file, meta_pageId, true);

        //build the tree
        FileScan fscan(relationName, bufMgrIn);
        try{
            RecordId scan_rid;
            while(1){
                fscan.scanNext(scan_rid);
                std::string recordStr = fscan.getRecord();
                const char* record = recordStr.c_str();
                if(attributeType == INTEGER){
                    int key = *((int*) (record + attrByteOffset));
                    std::cout << "Inserting : " << key << " into tree\n";
                    BTreeIndex::insertEntry(&key, scan_rid);
                }
                else if (attributeType == DOUBLE){
                    double key = *((double*) (record + attrByteOffset));
                    std::cout << "Inserting : " << key <<" into tree\n";
                    BTreeIndex::insertEntry(&key, scan_rid);
                }
                else{
                    char keyHolder[10];
                    char* src = (char*)(record + attrByteOffset);
                    strncpy(keyHolder, src, sizeof(keyHolder));
                    std::string key = std::string(keyHolder);
                    std::cout << "Inserting : " << key << " into tree\n";
                    BTreeIndex::insertEntry(&key, scan_rid);
                }
            }
        }
        catch(EndOfFileException e){
            scanExecuting = false;
            std::cout << "All records had been read in constructors\n";
        }

        std::cout << "File created: " << (*metapage).relationName << "\n";
        bufMgr->flushFile(file);
    }
    
    catch(FileExistsException ex){
        //check for metafile
        BlobFile* newfile = new BlobFile(outIndexName, false);
        file = (File*) newfile;
        IndexMetaInfo* metapage_info;
        Page *meta_page;
        bufMgrIn->readPage(file, headerPageNum, meta_page);
        metapage_info = (IndexMetaInfo*) meta_page;
        rootPageNum = (*metapage_info).rootPageNo;
        std::cout << (*metapage_info).relationName << "  "<< relationName<< "\n";
        if(relationName.compare((*metapage_info).relationName)==0 
            && (*metapage_info).attrByteOffset==attrByteOffset 
            && (*metapage_info).attrType==attrType){
            rootPageNum = (*metapage_info).rootPageNo;

            return;
        }
        else{
            throw BadIndexInfoException("Metafile does not match!!");
        }
        bufMgr->unPinPage(file, headerPageNum, false);

    }catch(FileNotFoundException ex){
        std::cout<<"You should not get FileNotFoundException in constructor!!\n";
    }	
}


// -----------------------------------------------------------------------------
// BTreeIndex::~BTreeIndex -- destructor
// -----------------------------------------------------------------------------
/**
	 * End any initialized scan, flush index file, after unpinning any pinned pages, from the buffer manager
	 * and delete file instance thereby closing the index file.
	 * Destructor should not throw any exceptions. All exceptions should be caught in here itself.
	 * */

BTreeIndex::~BTreeIndex()
{
    scanExecuting = false;

    try {
        bufMgr->unPinPage(file, currentPageNum, false);
    } catch (PageNotPinnedException e) {
    } catch (HashNotFoundException e) {
    }

    bufMgr->flushFile(file);
    file->~File();
}

// -----------------------------------------------------------------------------
// BTreeIndex::insertEntry
// -----------------------------------------------------------------------------

const void BTreeIndex::insertEntry(const void *key, const RecordId rid) 
{
    if(attributeType == INTEGER) {
        RIDKeyPair<int> rInt;
        rInt.set(rid, *((int *) key));
        int newValueInt;
        PageId newPidInt = 0;

        //call insertInto function recursively
        insertInto<int, LeafNodeInt, NonLeafNodeInt>
        (rInt, rootPageNum, 0, INTARRAYLEAFSIZE, INTARRAYNONLEAFSIZE, newValueInt, newPidInt);
        
        //if root got split
        if (newPidInt != 0)
            reArrangeRoot<int, NonLeafNodeInt>
            (newValueInt, newPidInt, INTARRAYNONLEAFSIZE);
            
    } else if (attributeType == DOUBLE) {
        RIDKeyPair<double> rDouble;
        rDouble.set(rid, *((double *) key));
        double newValueDouble;
        PageId newPidDouble = 0;

        //call recursive function
        insertInto<double, LeafNodeDouble, NonLeafNodeDouble>
                (rDouble, rootPageNum, 0, DOUBLEARRAYLEAFSIZE, DOUBLEARRAYNONLEAFSIZE, newValueDouble, newPidDouble);

        //if root got split
        if (newPidDouble != 0)
            reArrangeRoot<double, NonLeafNodeDouble>(newValueDouble, newPidDouble, DOUBLEARRAYNONLEAFSIZE);
    } else if (attributeType == STRING){
        RIDKeyPair<char[STRINGSIZE] > rString;
        rString.rid = rid;
        strncpy(rString.key, (char*) key, STRINGSIZE);
        char newValue[STRINGSIZE];
        PageId newPageId = 0;

        //call recursive function
        insertInto<char[STRINGSIZE], LeafNodeString, NonLeafNodeString>
                (rString, rootPageNum, 0, STRINGARRAYLEAFSIZE, STRINGARRAYNONLEAFSIZE, newValue, newPageId);

        //if root got split
        if (newPageId != 0)
            reArrangeRoot<char[STRINGSIZE], NonLeafNodeString>(newValue, newPageId, STRINGARRAYNONLEAFSIZE);
    }

}

/*
 * This function will search a position for a new entry.
 * If current page is a non-leaf page, it will find a position and call insertEntryRecursive.
 * If current page is leaf page, it will try to insert new entry.
 * */
template <class T, class T1, class T2>
void BTreeIndex::insertInto(RIDKeyPair<T> ridKeyPair,
                                      PageId pageId,
                                      bool isLeaf,
                                      int LEAFARRAYMAX,
                                      int NONLEAFARRAYMAX,
                                      T & newValue,
                                      PageId& newPageId)
{
    if(isLeaf){
        // Read the page
        Page* page;
        bufMgr->readPage(file, pageId, page);
        T1 *leafNode = (T1 *) page;

        // Find position
        int pos = 0;
        while(biggerThan<T> (ridKeyPair.key, leafNode->keyArray[pos])
              && leafNode->ridArray[pos].page_number != 0
              && pos < LEAFARRAYMAX)
            pos++;

        // Find last entry
        int last;
        for (last =0; last < LEAFARRAYMAX; last++)
            if(leafNode->ridArray[last].page_number == 0)
                break;


        if(last < LEAFARRAYMAX){
            // Not full
            for(int i=last; i>pos; i--){
                copy<T> (leafNode->keyArray[i], leafNode->keyArray[i - 1]);
                leafNode->ridArray[i] = leafNode->ridArray[i - 1];
            }
            copy<T> (leafNode->keyArray[pos], ridKeyPair.key);
            leafNode->ridArray[pos] = ridKeyPair.rid;
        } else {
            // Full, call split helper.
            leafSplitHelper<T,T1>(pos, last, LEAFARRAYMAX,
                    NONLEAFARRAYMAX, ridKeyPair,leafNode,newPageId,newValue);
        }

        leafOccupancy++;
        //unpin
        bufMgr->unPinPage(file, pageId, true);

    } else {
        // Non leaf
        // Read page
        Page* page;
        bufMgr->readPage(file, pageId, page);
        T2* nonLeafNode = (T2*) page;

        // Find pageArray position.
        int pos = 0;
        
        while(!smallerThan<T> (ridKeyPair.key, nonLeafNode->keyArray[pos])
              && nonLeafNode->pageNoArray[pos + 1] != 0
              && pos < NONLEAFARRAYMAX)
            pos++;

        // Index file is empty.
        if(nonLeafNode->pageNoArray[pos] == 0){
            PageId newPageId;
            Page* page;
            bufMgr->allocPage(file, newPageId, page);

            T1* leafNode = (T1*) page;
            for(int i=0; i < LEAFARRAYMAX; i++)
                leafNode->ridArray[i].page_number = 0;

            copy<T> (leafNode->keyArray[0], ridKeyPair.key);
            leafNode->ridArray[0] = ridKeyPair.rid;
            leafNode->rightSibPageNo = 0;
            nonLeafNode->pageNoArray[0] = newPageId;

            //unpin page
            bufMgr->unPinPage(file, newPageId, true);
            bufMgr->unPinPage(file, pageId, true);
            
            return;
        }

        // Call recursive function.
        bufMgr->unPinPage(file, pageId, false);
        T newChildValue;
        PageId newChildPageId = 0;
        insertInto<T, T1, T2>(ridKeyPair, nonLeafNode->pageNoArray[pos], nonLeafNode->level == 1,
                                        LEAFARRAYMAX, NONLEAFARRAYMAX, newChildValue, newChildPageId);

        // Check if child split.
        if(newChildPageId != 0){
            // If child split.
            bufMgr->readPage(file, pageId, page);
            T2* nonLeafNode = (T2*) page;

            // Find last entry.
            int last;
            for (last = 0; last < NONLEAFARRAYMAX; last++)
                if(nonLeafNode->pageNoArray[last + 1] == 0)
                    break;

            // Check if full.
            if(last < NONLEAFARRAYMAX){
                // Not full, just insert.
                for(int i=last; i>pos; i--){
                    copy<T> (nonLeafNode->keyArray[i], nonLeafNode->keyArray[i - 1]);
                    nonLeafNode->pageNoArray[i + 1] = nonLeafNode->pageNoArray[i];
                }
                copy<T> (nonLeafNode->keyArray[pos], newChildValue);
                nonLeafNode->pageNoArray[pos + 1] = newChildPageId;
            }
            else{
                // Full. Call split helper.
                nonLeafSplitHelper<T, T2>(pos, NONLEAFARRAYMAX, nonLeafNode,newPageId,newValue,
                                          newChildValue,newChildPageId);
            }

            nodeOccupancy++;
            bufMgr->unPinPage(file, pageId, true);
        }
    }
}

/*
 * This function is called when root node got split.
 * We need to create a new root page and link it to the old root page and new page.
 * */
template<class T, class T1>
void BTreeIndex::reArrangeRoot(T& newValue, PageId newPageId, int ARRAYMAX){
    PageId newRootPageId;
    Page *newRootPage;
    bufMgr->allocPage(file, newRootPageId, newRootPage);

    T1 *newRootNonLeafNode = (T1 *) newRootPage;
    for(int i=0;i<ARRAYMAX + 1; i++) newRootNonLeafNode->pageNoArray[i] = 0;
    copy<T> (newRootNonLeafNode->keyArray[0], newValue);
    newRootNonLeafNode->pageNoArray[0] = rootPageNum;
    newRootNonLeafNode->pageNoArray[1] = newPageId;
    newRootNonLeafNode->level = 0;
    rootPageNum = newRootPageId;
    bufMgr->unPinPage(file, newRootPageId, true);
}

/*
 * this helper help the leaf to split and return newPageId and new Value to parent
 * */
template <class T, class T1>
void BTreeIndex::leafSplitHelper(int pos, int last, int LEAFARRAYMAX,
                                 int NONLEAFARRAYMAX,
                                 RIDKeyPair<T> ridKeyPair,
                                 T1* leafNode,
                              	 PageId& newPageId,
                                 T & newValue)
{
    //full
    Page* newPage;
    bufMgr->allocPage(file, newPageId, newPage);
    T1 *newLeafNode = (T1 *) newPage;

    //tmp array
    T tmpKeyArray[LEAFARRAYMAX + 1];
    RecordId tmpRidArray[LEAFARRAYMAX + 1];

    //copy all new records to tmp
    for(int i=0; i < LEAFARRAYMAX; i++) {
        tmpRidArray[i] = leafNode->ridArray[i];
        copy<T> (tmpKeyArray[i], leafNode->keyArray[i]);
        leafNode->ridArray[i].page_number = 0;
    }

    for(int i= LEAFARRAYMAX; i > pos; i--){
        copy<T> (tmpKeyArray[i], tmpKeyArray[i-1]);
        tmpRidArray[i] = tmpRidArray[i-1];
    }
    copy<T> (tmpKeyArray[pos], ridKeyPair.key);
    tmpRidArray[pos] = ridKeyPair.rid;

    //reset new and old page
    for(int i=0; i < LEAFARRAYMAX; i++) {
        leafNode->ridArray[i].page_number = 0;
        newLeafNode->ridArray[i].page_number = 0;
    }

    //copy back
    for(int i=0; i< (LEAFARRAYMAX + 1) / 2; i++ ){
        copy<T> (leafNode->keyArray[i], tmpKeyArray[i]);
        leafNode->ridArray[i] = tmpRidArray[i];
    }
    for(int i= (LEAFARRAYMAX + 1) / 2; i < LEAFARRAYMAX + 1; i++){
        copy<T> (newLeafNode->keyArray[i - (LEAFARRAYMAX + 1) / 2], tmpKeyArray[i]);
        newLeafNode->ridArray[i - (LEAFARRAYMAX + 1) / 2] = tmpRidArray[i];
    }

    //link leaf node
    newLeafNode->rightSibPageNo = leafNode->rightSibPageNo;
    leafNode->rightSibPageNo = newPageId;

    //push up
    copy<T> (newValue, newLeafNode->keyArray[0]);

    //unpin
    bufMgr->unPinPage(file, newPageId, true);
};



/*
 * When the non Leaf need to split it call this nonLeafSplit Helper
 * */
template <class T, class T2>
void  BTreeIndex::nonLeafSplitHelper(int pos,
                                     int NONLEAFARRAYMAX,
                                     T2* nonLeafNode,
                                     PageId& newPageId,
                                     T & newValue,
                                     T&  newChildValue,
                                     PageId newChildPageId){
    //full, need split
    Page* newPage;
    bufMgr->allocPage(file, newPageId, newPage);
    T2* newNonLeafNode = (T2*) newPage;

    //tmp array
    T tmpKeyArray[NONLEAFARRAYMAX + 1];
    PageId tmpPageIdArray[NONLEAFARRAYMAX + 2];

    //copy to tmp
    for(int i=0; i < NONLEAFARRAYMAX; i++) {
        tmpPageIdArray[i] = nonLeafNode->pageNoArray[i];
        copy<T> (tmpKeyArray[i], nonLeafNode->keyArray[i]);
    }
    tmpPageIdArray[NONLEAFARRAYMAX + 1] = nonLeafNode->pageNoArray[NONLEAFARRAYMAX + 1];

    for(int i= NONLEAFARRAYMAX; i > pos; i--){
        copy<T> (tmpKeyArray[i], tmpKeyArray[i-1]);
        tmpPageIdArray[i+1] = tmpPageIdArray[i];
    }
    copy<T> (tmpKeyArray[pos], newChildValue);
    tmpPageIdArray[pos+1] = newChildPageId;

    //clear old and new page
    for(int i=0; i < NONLEAFARRAYMAX + 1; i++){
        nonLeafNode->pageNoArray[i] = 0;
        newNonLeafNode->pageNoArray[i] = 0;
    }

    //copy back
    for(int i=0; i< (NONLEAFARRAYMAX + 1) / 2; i++ ){
        copy<T> (nonLeafNode->keyArray[i], tmpKeyArray[i]);
        nonLeafNode->pageNoArray[i] = tmpPageIdArray[i];
    }
    nonLeafNode->pageNoArray[(NONLEAFARRAYMAX + 1) / 2] = tmpPageIdArray[(NONLEAFARRAYMAX + 1) / 2];

    for(int i= (NONLEAFARRAYMAX + 1) / 2 + 1; i < NONLEAFARRAYMAX + 1; i++){
        copy<T> (newNonLeafNode->keyArray[i - (NONLEAFARRAYMAX + 1) / 2 - 1], tmpKeyArray[i]);
        newNonLeafNode->pageNoArray[i - (NONLEAFARRAYMAX + 1) / 2 - 1 ] = tmpPageIdArray[i];
    }
    newNonLeafNode->pageNoArray[NONLEAFARRAYMAX + 1 - (NONLEAFARRAYMAX + 1) / 2 - 1 ] = tmpPageIdArray[
            NONLEAFARRAYMAX + 1];

    //level
    newNonLeafNode->level = nonLeafNode->level;

    //push up
    copy<T> (newValue, tmpKeyArray[(NONLEAFARRAYMAX + 1) / 2]);

    //unpin
    bufMgr->unPinPage(file, newPageId, true);
};

// -----------------------------------------------------------------------------
// BTreeIndex::startScan
// -----------------------------------------------------------------------------


template<class T, class T1>
void BTreeIndex::startScanHelper(T lowValParm,
                                                                 T highValParm,
                                                                 int NONLEAFARRAYMAX)
{
        if(biggerThan<T> (lowValParm,  highValParm))
                throw BadScanrangeException();

        //find first one
        currentPageNum = rootPageNum;
        bufMgr->readPage(file, currentPageNum, currentPageData);
        T1* nonLeafNode = (T1*) currentPageData;

        int pos = 0;
        while(nonLeafNode->level != 1) {
                // If current level is not 1, then next page is not leaf page.
                // Still need to go to next level.
                pos = 0;
                while(!smallerThan<T> (lowValParm, nonLeafNode->keyArray[pos])
                            && nonLeafNode->pageNoArray[pos + 1] != 0
                            && pos < NONLEAFARRAYMAX)
                        pos++;
                PageId nextPageId = nonLeafNode->pageNoArray[pos];
                bufMgr->readPage(file, nextPageId, currentPageData);
                bufMgr->unPinPage(file, currentPageNum, false);
                currentPageNum = nextPageId;
                nonLeafNode = (T1*) currentPageData;
        }

        // This page is level 1, which means next page is leaf node.
        pos = 0;
        while(!smallerThan<T> (lowValParm, nonLeafNode->keyArray[pos])
                    && nonLeafNode->pageNoArray[pos + 1] != 0
                    && pos < NONLEAFARRAYMAX)
                pos++;
        PageId nextPageId = nonLeafNode->pageNoArray[pos];
        bufMgr->readPage(file, nextPageId, currentPageData);
        bufMgr->unPinPage(file, currentPageNum, false);
        currentPageNum = nextPageId;
        nextEntry = 0;
}

const void BTreeIndex::startScan(const void* lowValParm,
				   const Operator lowOpParm,
				   const void* highValParm,
				   const Operator highOpParm)
{
    scanExecuting = true;
    if (lowOpParm != GT && lowOpParm !=GTE ){
        throw BadOpcodesException();
    }
    if(highOpParm != LT && highOpParm != LTE){
        throw BadOpcodesException();
    }

    this->lowOp = lowOpParm;
    this->highOp = highOpParm;

    // Call helper to do the work.
    if(attributeType == INTEGER) {
        this->lowValInt = *((int*) lowValParm);
        this->highValInt = *((int *) highValParm);
        startScanHelper<int, NonLeafNodeInt>(*((int*) lowValParm), *((int *) highValParm), INTARRAYNONLEAFSIZE);
    } else if (attributeType == DOUBLE) {
        this->lowValDouble = *((double*) lowValParm);
        this->highValDouble = *((double *) highValParm);
        startScanHelper<double , NonLeafNodeDouble>(*((double *) lowValParm), *((double *) highValParm), DOUBLEARRAYNONLEAFSIZE);
    } else {
        strncpy((char*) this->lowValString.c_str(), (char *)lowValParm, STRINGSIZE-1);
        strncpy(lowValChar, (char *)lowValParm, STRINGSIZE);
        strncpy((char*) this->highValString.c_str(), (char *)highValParm, STRINGSIZE-1);
        strncpy(highValChar, (char *)highValParm, STRINGSIZE);
        startScanHelper<char[STRINGSIZE], NonLeafNodeString> (lowValChar, highValChar, STRINGARRAYNONLEAFSIZE);
    }

}


// -----------------------------------------------------------------------------
// BTreeIndex::scanNext
// -----------------------------------------------------------------------------

template <class T, class T1>
void BTreeIndex::scanNextHelper(RecordId &outRid, T lowVal, T highVal, int ARRAYMAX)
{
    T1* leafNode;
    while(1){
        leafNode = (T1*) currentPageData;

        // Go to next page.
        if(leafNode->ridArray[nextEntry].page_number == 0 || nextEntry == ARRAYMAX) {
            PageId nextPageNum = leafNode->rightSibPageNo;
            if(nextPageNum == 0){
                // Next page is 0, scan finish.
                bufMgr->unPinPage(file, currentPageNum, false);
                throw IndexScanCompletedException();
            }

            bufMgr->unPinPage(file, currentPageNum, false);
            currentPageNum = nextPageNum;

            bufMgr->readPage(file, currentPageNum, currentPageData);
            nextEntry = 0;
            continue;
        }

        // Do not satisfy.
        if((lowOp==GT && !biggerThan<T> (leafNode->keyArray[nextEntry], lowVal) )
           || (lowOp==GTE && smallerThan<T> (leafNode->keyArray[nextEntry], lowVal))
           ) {
            nextEntry++;
            continue;
        }

        // Value bigger that high value, scan end.
        if((highOp==LT && !smallerThan<T> (leafNode->keyArray[nextEntry],  highVal))
           || (highOp==LTE && biggerThan<T> (leafNode->keyArray[nextEntry], highVal) ))
            throw IndexScanCompletedException();

        // Got a record.
        outRid = leafNode->ridArray[nextEntry];
        nextEntry++;
        return ;
    }
}

const void BTreeIndex::scanNext(RecordId& outRid) 
{
    if(!scanExecuting)
        throw ScanNotInitializedException();

    if(attributeType == INTEGER)
        scanNextHelper<int, LeafNodeInt> (outRid, lowValInt, highValInt, INTARRAYLEAFSIZE);
    else if (attributeType == DOUBLE)
        scanNextHelper<double, LeafNodeDouble> (outRid, lowValDouble, highValDouble, DOUBLEARRAYLEAFSIZE);
    else
        scanNextHelper<char[STRINGSIZE], LeafNodeString> (outRid, lowValChar, highValChar, STRINGARRAYLEAFSIZE);
}

// -----------------------------------------------------------------------------
// BTreeIndex::endScan
// -----------------------------------------------------------------------------
//
const void BTreeIndex::endScan() 
{
    if(!scanExecuting)
        throw ScanNotInitializedException();

    try {
        bufMgr->unPinPage(file, currentPageNum, false);
    } catch (PageNotPinnedException e) {
    } catch (HashNotFoundException e) {
    }
    
    scanExecuting = false;
}

}
