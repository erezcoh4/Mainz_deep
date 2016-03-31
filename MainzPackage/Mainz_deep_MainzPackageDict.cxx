// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME Mainz_deep_MainzPackageDict

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
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

// Since CINT ignores the std namespace, we need to do so in this file.
namespace std {} using namespace std;

// Header files passed as explicit arguments
#include "CK_RungeKutta.h"
#include "DipolB.h"
#include "FringeFall.h"
#include "Magnets.h"
#include "Matrix3D.h"
#include "TCalcDeep.h"
#include "Tqspin.h"
#include "TqspinSpectrometer.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static TClass *TqspinSpectrometer_Dictionary();
   static void TqspinSpectrometer_TClassManip(TClass*);
   static void delete_TqspinSpectrometer(void *p);
   static void deleteArray_TqspinSpectrometer(void *p);
   static void destruct_TqspinSpectrometer(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TqspinSpectrometer*)
   {
      ::TqspinSpectrometer *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::TqspinSpectrometer));
      static ::ROOT::TGenericClassInfo 
         instance("TqspinSpectrometer", "TqspinSpectrometer.h", 44,
                  typeid(::TqspinSpectrometer), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &TqspinSpectrometer_Dictionary, isa_proxy, 4,
                  sizeof(::TqspinSpectrometer) );
      instance.SetDelete(&delete_TqspinSpectrometer);
      instance.SetDeleteArray(&deleteArray_TqspinSpectrometer);
      instance.SetDestructor(&destruct_TqspinSpectrometer);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TqspinSpectrometer*)
   {
      return GenerateInitInstanceLocal((::TqspinSpectrometer*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::TqspinSpectrometer*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *TqspinSpectrometer_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::TqspinSpectrometer*)0x0)->GetClass();
      TqspinSpectrometer_TClassManip(theClass);
   return theClass;
   }

   static void TqspinSpectrometer_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *CK_RungeKutta_Dictionary();
   static void CK_RungeKutta_TClassManip(TClass*);
   static void *new_CK_RungeKutta(void *p = 0);
   static void *newArray_CK_RungeKutta(Long_t size, void *p);
   static void delete_CK_RungeKutta(void *p);
   static void deleteArray_CK_RungeKutta(void *p);
   static void destruct_CK_RungeKutta(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::CK_RungeKutta*)
   {
      ::CK_RungeKutta *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::CK_RungeKutta));
      static ::ROOT::TGenericClassInfo 
         instance("CK_RungeKutta", "CK_RungeKutta.h", 30,
                  typeid(::CK_RungeKutta), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &CK_RungeKutta_Dictionary, isa_proxy, 4,
                  sizeof(::CK_RungeKutta) );
      instance.SetNew(&new_CK_RungeKutta);
      instance.SetNewArray(&newArray_CK_RungeKutta);
      instance.SetDelete(&delete_CK_RungeKutta);
      instance.SetDeleteArray(&deleteArray_CK_RungeKutta);
      instance.SetDestructor(&destruct_CK_RungeKutta);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::CK_RungeKutta*)
   {
      return GenerateInitInstanceLocal((::CK_RungeKutta*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::CK_RungeKutta*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *CK_RungeKutta_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::CK_RungeKutta*)0x0)->GetClass();
      CK_RungeKutta_TClassManip(theClass);
   return theClass;
   }

   static void CK_RungeKutta_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *Tqspin_Dictionary();
   static void Tqspin_TClassManip(TClass*);
   static void *new_Tqspin(void *p = 0);
   static void *newArray_Tqspin(Long_t size, void *p);
   static void delete_Tqspin(void *p);
   static void deleteArray_Tqspin(void *p);
   static void destruct_Tqspin(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::Tqspin*)
   {
      ::Tqspin *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::Tqspin));
      static ::ROOT::TGenericClassInfo 
         instance("Tqspin", "Tqspin.h", 39,
                  typeid(::Tqspin), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &Tqspin_Dictionary, isa_proxy, 4,
                  sizeof(::Tqspin) );
      instance.SetNew(&new_Tqspin);
      instance.SetNewArray(&newArray_Tqspin);
      instance.SetDelete(&delete_Tqspin);
      instance.SetDeleteArray(&deleteArray_Tqspin);
      instance.SetDestructor(&destruct_Tqspin);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::Tqspin*)
   {
      return GenerateInitInstanceLocal((::Tqspin*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::Tqspin*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *Tqspin_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::Tqspin*)0x0)->GetClass();
      Tqspin_TClassManip(theClass);
   return theClass;
   }

   static void Tqspin_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *Dipol_Dictionary();
   static void Dipol_TClassManip(TClass*);
   static void *new_Dipol(void *p = 0);
   static void *newArray_Dipol(Long_t size, void *p);
   static void delete_Dipol(void *p);
   static void deleteArray_Dipol(void *p);
   static void destruct_Dipol(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::Dipol*)
   {
      ::Dipol *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::Dipol));
      static ::ROOT::TGenericClassInfo 
         instance("Dipol", "Magnets.h", 97,
                  typeid(::Dipol), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &Dipol_Dictionary, isa_proxy, 4,
                  sizeof(::Dipol) );
      instance.SetNew(&new_Dipol);
      instance.SetNewArray(&newArray_Dipol);
      instance.SetDelete(&delete_Dipol);
      instance.SetDeleteArray(&deleteArray_Dipol);
      instance.SetDestructor(&destruct_Dipol);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::Dipol*)
   {
      return GenerateInitInstanceLocal((::Dipol*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::Dipol*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *Dipol_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::Dipol*)0x0)->GetClass();
      Dipol_TClassManip(theClass);
   return theClass;
   }

   static void Dipol_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *TCalcDeep_Dictionary();
   static void TCalcDeep_TClassManip(TClass*);
   static void *new_TCalcDeep(void *p = 0);
   static void *newArray_TCalcDeep(Long_t size, void *p);
   static void delete_TCalcDeep(void *p);
   static void deleteArray_TCalcDeep(void *p);
   static void destruct_TCalcDeep(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TCalcDeep*)
   {
      ::TCalcDeep *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::TCalcDeep));
      static ::ROOT::TGenericClassInfo 
         instance("TCalcDeep", "TCalcDeep.h", 28,
                  typeid(::TCalcDeep), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &TCalcDeep_Dictionary, isa_proxy, 4,
                  sizeof(::TCalcDeep) );
      instance.SetNew(&new_TCalcDeep);
      instance.SetNewArray(&newArray_TCalcDeep);
      instance.SetDelete(&delete_TCalcDeep);
      instance.SetDeleteArray(&deleteArray_TCalcDeep);
      instance.SetDestructor(&destruct_TCalcDeep);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TCalcDeep*)
   {
      return GenerateInitInstanceLocal((::TCalcDeep*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::TCalcDeep*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *TCalcDeep_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::TCalcDeep*)0x0)->GetClass();
      TCalcDeep_TClassManip(theClass);
   return theClass;
   }

   static void TCalcDeep_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrapper around operator delete
   static void delete_TqspinSpectrometer(void *p) {
      delete ((::TqspinSpectrometer*)p);
   }
   static void deleteArray_TqspinSpectrometer(void *p) {
      delete [] ((::TqspinSpectrometer*)p);
   }
   static void destruct_TqspinSpectrometer(void *p) {
      typedef ::TqspinSpectrometer current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::TqspinSpectrometer

namespace ROOT {
   // Wrappers around operator new
   static void *new_CK_RungeKutta(void *p) {
      return  p ? new(p) ::CK_RungeKutta : new ::CK_RungeKutta;
   }
   static void *newArray_CK_RungeKutta(Long_t nElements, void *p) {
      return p ? new(p) ::CK_RungeKutta[nElements] : new ::CK_RungeKutta[nElements];
   }
   // Wrapper around operator delete
   static void delete_CK_RungeKutta(void *p) {
      delete ((::CK_RungeKutta*)p);
   }
   static void deleteArray_CK_RungeKutta(void *p) {
      delete [] ((::CK_RungeKutta*)p);
   }
   static void destruct_CK_RungeKutta(void *p) {
      typedef ::CK_RungeKutta current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::CK_RungeKutta

namespace ROOT {
   // Wrappers around operator new
   static void *new_Tqspin(void *p) {
      return  p ? new(p) ::Tqspin : new ::Tqspin;
   }
   static void *newArray_Tqspin(Long_t nElements, void *p) {
      return p ? new(p) ::Tqspin[nElements] : new ::Tqspin[nElements];
   }
   // Wrapper around operator delete
   static void delete_Tqspin(void *p) {
      delete ((::Tqspin*)p);
   }
   static void deleteArray_Tqspin(void *p) {
      delete [] ((::Tqspin*)p);
   }
   static void destruct_Tqspin(void *p) {
      typedef ::Tqspin current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::Tqspin

namespace ROOT {
   // Wrappers around operator new
   static void *new_Dipol(void *p) {
      return  p ? new(p) ::Dipol : new ::Dipol;
   }
   static void *newArray_Dipol(Long_t nElements, void *p) {
      return p ? new(p) ::Dipol[nElements] : new ::Dipol[nElements];
   }
   // Wrapper around operator delete
   static void delete_Dipol(void *p) {
      delete ((::Dipol*)p);
   }
   static void deleteArray_Dipol(void *p) {
      delete [] ((::Dipol*)p);
   }
   static void destruct_Dipol(void *p) {
      typedef ::Dipol current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::Dipol

namespace ROOT {
   // Wrappers around operator new
   static void *new_TCalcDeep(void *p) {
      return  p ? new(p) ::TCalcDeep : new ::TCalcDeep;
   }
   static void *newArray_TCalcDeep(Long_t nElements, void *p) {
      return p ? new(p) ::TCalcDeep[nElements] : new ::TCalcDeep[nElements];
   }
   // Wrapper around operator delete
   static void delete_TCalcDeep(void *p) {
      delete ((::TCalcDeep*)p);
   }
   static void deleteArray_TCalcDeep(void *p) {
      delete [] ((::TCalcDeep*)p);
   }
   static void destruct_TCalcDeep(void *p) {
      typedef ::TCalcDeep current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::TCalcDeep

namespace {
  void TriggerDictionaryInitialization_libMainz_deep_MainzPackage_Impl() {
    static const char* headers[] = {
"CK_RungeKutta.h",
"DipolB.h",
"FringeFall.h",
"Magnets.h",
"Matrix3D.h",
"TCalcDeep.h",
"Tqspin.h",
"TqspinSpectrometer.h",
0
    };
    static const char* includePaths[] = {
"/Users/erezcohen/larlite/UserDev/mySoftware",
"/usr/local/Cellar/root6/6.06.02/include/root",
"/Users/erezcohen/larlite/UserDev/Mainz_deep/MainzPackage/",
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "libMainz_deep_MainzPackage dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
class __attribute__((annotate("$clingAutoload$TCalcDeep.h")))  TqspinSpectrometer;
class __attribute__((annotate("$clingAutoload$CK_RungeKutta.h")))  CK_RungeKutta;
class __attribute__((annotate("$clingAutoload$TCalcDeep.h")))  Tqspin;
class __attribute__((annotate("$clingAutoload$DipolB.h")))  Dipol;
class __attribute__((annotate("$clingAutoload$TCalcDeep.h")))  TCalcDeep;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "libMainz_deep_MainzPackage dictionary payload"

#ifndef G__VECTOR_HAS_CLASS_ITERATOR
  #define G__VECTOR_HAS_CLASS_ITERATOR 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
#include "CK_RungeKutta.h"
#include "DipolB.h"
#include "FringeFall.h"
#include "Magnets.h"
#include "Matrix3D.h"
#include "TCalcDeep.h"
#include "Tqspin.h"
#include "TqspinSpectrometer.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
"CK_RungeKutta", payloadCode, "@",
"Dipol", payloadCode, "@",
"TCalcDeep", payloadCode, "@",
"Tqspin", payloadCode, "@",
"TqspinSpectrometer", payloadCode, "@",
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("libMainz_deep_MainzPackage",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_libMainz_deep_MainzPackage_Impl, {}, classesHeaders);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_libMainz_deep_MainzPackage_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_libMainz_deep_MainzPackage() {
  TriggerDictionaryInitialization_libMainz_deep_MainzPackage_Impl();
}
