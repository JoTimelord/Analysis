// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME selectordict
#define R__NO_DEPRECATION

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "RConfig.h"
#include "TClass.h"
#include "TDictAttributeMap.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TInterpreter.h"
#include "TVirtualMutex.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TIsAProxy.h"
#include "TFileMergeInfo.h"
#include <algorithm>
#include "TCollectionProxyInfo.h"
/*******************************************************************/

#include "TDataMember.h"

// The generated code does not explicitly qualifies STL entities
namespace std {} using namespace std;

// Header files passed as explicit arguments
#include "Selector.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static void *new_Selector(void *p = 0);
   static void *newArray_Selector(Long_t size, void *p);
   static void delete_Selector(void *p);
   static void deleteArray_Selector(void *p);
   static void destruct_Selector(void *p);
   static void streamer_Selector(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::Selector*)
   {
      ::Selector *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::Selector >(0);
      static ::ROOT::TGenericClassInfo 
         instance("Selector", ::Selector::Class_Version(), "Selector.h", 22,
                  typeid(::Selector), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::Selector::Dictionary, isa_proxy, 16,
                  sizeof(::Selector) );
      instance.SetNew(&new_Selector);
      instance.SetNewArray(&newArray_Selector);
      instance.SetDelete(&delete_Selector);
      instance.SetDeleteArray(&deleteArray_Selector);
      instance.SetDestructor(&destruct_Selector);
      instance.SetStreamerFunc(&streamer_Selector);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::Selector*)
   {
      return GenerateInitInstanceLocal((::Selector*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::Selector*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr Selector::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *Selector::Class_Name()
{
   return "Selector";
}

//______________________________________________________________________________
const char *Selector::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::Selector*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int Selector::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::Selector*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *Selector::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::Selector*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *Selector::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::Selector*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void Selector::Streamer(TBuffer &R__b)
{
   // Stream an object of class Selector.

   TSelector::Streamer(R__b);
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_Selector(void *p) {
      return  p ? new(p) ::Selector : new ::Selector;
   }
   static void *newArray_Selector(Long_t nElements, void *p) {
      return p ? new(p) ::Selector[nElements] : new ::Selector[nElements];
   }
   // Wrapper around operator delete
   static void delete_Selector(void *p) {
      delete ((::Selector*)p);
   }
   static void deleteArray_Selector(void *p) {
      delete [] ((::Selector*)p);
   }
   static void destruct_Selector(void *p) {
      typedef ::Selector current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_Selector(TBuffer &buf, void *obj) {
      ((::Selector*)obj)->::Selector::Streamer(buf);
   }
} // end of namespace ROOT for class ::Selector

namespace {
  void TriggerDictionaryInitialization_selectordict_Impl() {
    static const char* headers[] = {
"Selector.h",
0
    };
    static const char* includePaths[] = {
"/cvmfs/cms.cern.ch/slc7_amd64_gcc900/lcg/root/6.22.08-8d9ab2b279c3f35e6100d909611c3c2f/include/",
"/home/users/joytzphysics/Analysis/tools/",
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "selectordict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_AutoLoading_Map;
class __attribute__((annotate("$clingAutoload$Selector.h")))  Selector;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "selectordict dictionary payload"


#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
#include "Selector.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[] = {
"Selector", payloadCode, "@",
nullptr
};
    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("selectordict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_selectordict_Impl, {}, classesHeaders, /*hasCxxModule*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_selectordict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_selectordict() {
  TriggerDictionaryInitialization_selectordict_Impl();
}
