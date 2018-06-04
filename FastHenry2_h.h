

/* this ALWAYS GENERATED file contains the definitions for the interfaces */


 /* File created by MIDL compiler version 7.00.0555 */
/* at Mon Jun 04 16:19:29 2018
 */
/* Compiler settings for FastHenry2.odl:
    Oicf, W1, Zp8, env=Win32 (32b run), target_arch=X86 7.00.0555 
    protocol : dce , ms_ext, c_ext, robust
    error checks: allocation ref bounds_check enum stub_data 
    VC __declspec() decoration level: 
         __declspec(uuid()), __declspec(selectany), __declspec(novtable)
         DECLSPEC_UUID(), MIDL_INTERFACE()
*/
/* @@MIDL_FILE_HEADING(  ) */

#pragma warning( disable: 4049 )  /* more than 64k source lines */


/* verify that the <rpcndr.h> version is high enough to compile this file*/
#ifndef __REQUIRED_RPCNDR_H_VERSION__
#define __REQUIRED_RPCNDR_H_VERSION__ 475
#endif

#include "rpc.h"
#include "rpcndr.h"

#ifndef __RPCNDR_H_VERSION__
#error this stub requires an updated version of <rpcndr.h>
#endif // __RPCNDR_H_VERSION__


#ifndef __FastHenry2_h_h__
#define __FastHenry2_h_h__

#if defined(_MSC_VER) && (_MSC_VER >= 1020)
#pragma once
#endif

/* Forward Declarations */ 

#ifndef __IFastHenry_FWD_DEFINED__
#define __IFastHenry_FWD_DEFINED__
typedef interface IFastHenry IFastHenry;
#endif 	/* __IFastHenry_FWD_DEFINED__ */


#ifndef __Document_FWD_DEFINED__
#define __Document_FWD_DEFINED__

#ifdef __cplusplus
typedef class Document Document;
#else
typedef struct Document Document;
#endif /* __cplusplus */

#endif 	/* __Document_FWD_DEFINED__ */


#ifdef __cplusplus
extern "C"{
#endif 



#ifndef __FastHenry_LIBRARY_DEFINED__
#define __FastHenry_LIBRARY_DEFINED__

/* library FastHenry */
/* [version][uuid] */ 


DEFINE_GUID(LIBID_FastHenry,0xF41B20B2,0x4CFB,0x41A9,0xB9,0xE9,0x9E,0x92,0x9D,0xE0,0x11,0xDF);

#ifndef __IFastHenry_DISPINTERFACE_DEFINED__
#define __IFastHenry_DISPINTERFACE_DEFINED__

/* dispinterface IFastHenry */
/* [uuid] */ 


DEFINE_GUID(DIID_IFastHenry,0x3BAD3C03,0x3868,0x4C3C,0xB2,0x65,0xBD,0xAE,0x2A,0x4F,0xA4,0xF5);

#if defined(__cplusplus) && !defined(CINTERFACE)

    MIDL_INTERFACE("3BAD3C03-3868-4C3C-B265-BDAE2A4FA4F5")
    IFastHenry : public IDispatch
    {
    };
    
#else 	/* C style interface */

    typedef struct IFastHenryVtbl
    {
        BEGIN_INTERFACE
        
        HRESULT ( STDMETHODCALLTYPE *QueryInterface )( 
            IFastHenry * This,
            /* [in] */ REFIID riid,
            /* [annotation][iid_is][out] */ 
            __RPC__deref_out  void **ppvObject);
        
        ULONG ( STDMETHODCALLTYPE *AddRef )( 
            IFastHenry * This);
        
        ULONG ( STDMETHODCALLTYPE *Release )( 
            IFastHenry * This);
        
        HRESULT ( STDMETHODCALLTYPE *GetTypeInfoCount )( 
            IFastHenry * This,
            /* [out] */ UINT *pctinfo);
        
        HRESULT ( STDMETHODCALLTYPE *GetTypeInfo )( 
            IFastHenry * This,
            /* [in] */ UINT iTInfo,
            /* [in] */ LCID lcid,
            /* [out] */ ITypeInfo **ppTInfo);
        
        HRESULT ( STDMETHODCALLTYPE *GetIDsOfNames )( 
            IFastHenry * This,
            /* [in] */ REFIID riid,
            /* [size_is][in] */ LPOLESTR *rgszNames,
            /* [range][in] */ UINT cNames,
            /* [in] */ LCID lcid,
            /* [size_is][out] */ DISPID *rgDispId);
        
        /* [local] */ HRESULT ( STDMETHODCALLTYPE *Invoke )( 
            IFastHenry * This,
            /* [in] */ DISPID dispIdMember,
            /* [in] */ REFIID riid,
            /* [in] */ LCID lcid,
            /* [in] */ WORD wFlags,
            /* [out][in] */ DISPPARAMS *pDispParams,
            /* [out] */ VARIANT *pVarResult,
            /* [out] */ EXCEPINFO *pExcepInfo,
            /* [out] */ UINT *puArgErr);
        
        END_INTERFACE
    } IFastHenryVtbl;

    interface IFastHenry
    {
        CONST_VTBL struct IFastHenryVtbl *lpVtbl;
    };

    

#ifdef COBJMACROS


#define IFastHenry_QueryInterface(This,riid,ppvObject)	\
    ( (This)->lpVtbl -> QueryInterface(This,riid,ppvObject) ) 

#define IFastHenry_AddRef(This)	\
    ( (This)->lpVtbl -> AddRef(This) ) 

#define IFastHenry_Release(This)	\
    ( (This)->lpVtbl -> Release(This) ) 


#define IFastHenry_GetTypeInfoCount(This,pctinfo)	\
    ( (This)->lpVtbl -> GetTypeInfoCount(This,pctinfo) ) 

#define IFastHenry_GetTypeInfo(This,iTInfo,lcid,ppTInfo)	\
    ( (This)->lpVtbl -> GetTypeInfo(This,iTInfo,lcid,ppTInfo) ) 

#define IFastHenry_GetIDsOfNames(This,riid,rgszNames,cNames,lcid,rgDispId)	\
    ( (This)->lpVtbl -> GetIDsOfNames(This,riid,rgszNames,cNames,lcid,rgDispId) ) 

#define IFastHenry_Invoke(This,dispIdMember,riid,lcid,wFlags,pDispParams,pVarResult,pExcepInfo,puArgErr)	\
    ( (This)->lpVtbl -> Invoke(This,dispIdMember,riid,lcid,wFlags,pDispParams,pVarResult,pExcepInfo,puArgErr) ) 

#endif /* COBJMACROS */


#endif 	/* C style interface */


#endif 	/* __IFastHenry_DISPINTERFACE_DEFINED__ */


DEFINE_GUID(CLSID_Document,0x511D52AD,0x9892,0x4DF2,0xB1,0x90,0xBD,0x7D,0xFD,0x8E,0x73,0x22);

#ifdef __cplusplus

class DECLSPEC_UUID("511D52AD-9892-4DF2-B190-BD7DFD8E7322")
Document;
#endif
#endif /* __FastHenry_LIBRARY_DEFINED__ */

/* Additional Prototypes for ALL interfaces */

/* end of Additional Prototypes */

#ifdef __cplusplus
}
#endif

#endif


