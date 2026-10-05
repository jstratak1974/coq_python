#!/usr/bin/env python3
"""Standalone Tkinter GUI for COQ detector characterization and QA.

Launch: python coq_python.py
Dependencies: python -m pip install -r requirements.txt
"""
from __future__ import annotations

import argparse
import csv
import json
import queue
import threading
import traceback
from dataclasses import asdict
from pathlib import Path

import numpy as np
import tkinter as tk
from tkinter import ttk, filedialog, messagebox
from tkinter.scrolledtext import ScrolledText

from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.widgets import RectangleSelector
from matplotlib.patches import Rectangle

import coq_core as core


class COQApp:
    def __init__(self, root):
        self.root=root
        root.title("COQ Python | Detector characterization and quality assurance")
        root.geometry("1480x940");root.minsize(1280,800)
        root.columnconfigure(0,weight=1);root.rowconfigure(1,weight=1)
        self.images=[];self.active=0;self.rois={};self.response=core.Response()
        self.mtf_result=None;self.nps_result=None;self.current_result=None
        self.events=queue.Queue();self.busy=False;self.controls=[];self.fields={}
        self.tables=core.AGDTables();self.result_context={};self.auto_scale=True
        self.status=tk.StringVar(value="Open original detector QC images, or select a synthetic demonstration.")
        self.modality=tk.StringVar(value="Radiology DR")
        self.linearize=tk.BooleanVar(value=False);self.rescale=tk.BooleanVar(value=False)
        self.roi_slot=tk.StringVar(value="Analysis");self.all_nps=tk.BooleanVar(value=True)
        self.chart_view=tk.StringVar(value="Curve")
        self._build()
        self.root.after(100,self._poll)
        self.root.protocol("WM_DELETE_WINDOW",self.close)

    def _button(self,parent,text,command,run=False,**kwargs):
        b=ttk.Button(parent,text=text,command=lambda:self.safe(command),**kwargs)
        if run:self.controls.append(b)
        return b

    def safe(self,fun):
        try:fun()
        except Exception as e:
            self.status.set(str(e));messagebox.showerror("COQ Python",str(e),parent=self.root)

    def _field(self,parent,label,key,default="",choices=None,width=15):
        row=ttk.Frame(parent);row.pack(fill="x",pady=3)
        ttk.Label(row,text=label,width=31).pack(side="left")
        var=tk.StringVar(value=str(default));self.fields[key]=var
        widget=ttk.Combobox(row,textvariable=var,values=choices,state="readonly",width=width) if choices else ttk.Entry(row,textvariable=var,width=width)
        widget.pack(side="left",fill="x",expand=True)
        return widget

    def value(self,key,cast=float):return cast(self.fields[key].get())
    def optional(self,key):return None if not self.fields[key].get().strip() else self.value(key)

    def _note(self,parent,text):
        ttk.Label(parent,text=text,wraplength=480,justify="left",foreground="#526173").pack(fill="x",pady=(4,7))

    def _tab(self,title):
        outer=ttk.Frame(self.notebook);self.notebook.add(outer,text=title)
        canvas=tk.Canvas(outer,highlightthickness=0,width=510)
        scroll=ttk.Scrollbar(outer,orient="vertical",command=canvas.yview)
        canvas.configure(yscrollcommand=scroll.set);scroll.pack(side="right",fill="y");canvas.pack(side="left",fill="both",expand=True)
        frame=ttk.Frame(canvas,padding=10);item=canvas.create_window(0,0,window=frame,anchor="nw")
        frame.bind("<Configure>",lambda _:canvas.configure(scrollregion=canvas.bbox("all")))
        canvas.bind("<Configure>",lambda e:canvas.itemconfigure(item,width=e.width))
        return frame

    def _build(self):
        style=ttk.Style(self.root)
        if "clam" in style.theme_names():style.theme_use("clam")
        style.configure("TButton",padding=(7,5));style.configure("TNotebook.Tab",padding=(7,6))
        style.configure("Title.TLabel",font=("TkDefaultFont",15,"bold"))
        bar=ttk.Frame(self.root,padding=(12,8));bar.grid(row=0,column=0,sticky="ew")
        ttk.Label(bar,text="COQ Python",style="Title.TLabel").pack(side="left")
        ttk.Label(bar,text="  Physical characterization / detector QA",foreground="#526173").pack(side="left")
        ttk.Combobox(bar,textvariable=self.modality,state="readonly",width=20,
                     values=("Radiology DR","Radiology CR","Mammography DR","Mammography CR","Angiography DR")).pack(side="right")
        self._button(bar,"Help",self.help).pack(side="right",padx=6)
        paned=ttk.Panedwindow(self.root,orient="horizontal");paned.grid(row=1,column=0,sticky="nsew",padx=10)
        left=ttk.Frame(paned,width=550);right=ttk.Frame(paned)
        paned.add(left,weight=0);paned.add(right,weight=1)
        files=ttk.LabelFrame(left,text="Images and acquisition geometry",padding=6);files.pack(fill="x")
        tools=ttk.Frame(files);tools.pack(fill="x")
        self._button(tools,"Open images",self.open_images).pack(side="left")
        self._button(tools,"Import RAW",self.raw_dialog).pack(side="left",padx=4)
        self._button(tools,"Remove",self.remove_image).pack(side="left")
        self._button(tools,"Clear",self.clear_images).pack(side="left",padx=4)
        ttk.Checkbutton(files,text="Apply DICOM rescale slope/intercept when opening",variable=self.rescale).pack(anchor="w")
        listframe=ttk.Frame(files);listframe.pack(fill="x",pady=4)
        self.image_list=tk.Listbox(listframe,height=4,exportselection=False)
        self.image_list.pack(side="left",fill="x",expand=True)
        scroll=ttk.Scrollbar(listframe,orient="vertical",command=self.image_list.yview);scroll.pack(side="right",fill="y")
        self.image_list.configure(yscrollcommand=scroll.set);self.image_list.bind("<<ListboxSelect>>",self.choose_image)
        row=ttk.Frame(files);row.pack(fill="x")
        self._button(row,"Move up",lambda:self.move_image(-1)).pack(side="left")
        self._button(row,"Move down",lambda:self.move_image(1)).pack(side="left",padx=4)
        self._button(row,"Metadata",self.metadata).pack(side="left")
        for label,key in (("dy (mm)","dy"),("dx (mm)","dx")):
            ttk.Label(row,text=label).pack(side="left",padx=(5,2));self.fields[key]=tk.StringVar(value="")
            ttk.Entry(row,textvariable=self.fields[key],width=7).pack(side="left")
        demo=ttk.Frame(files);demo.pack(fill="x",pady=(5,0))
        self.demo_var=tk.StringVar(value="Slanted edge")
        ttk.Combobox(demo,textvariable=self.demo_var,width=22,state="readonly",
                     values=("Slanted edge","Slanted slit","Flat fields","Response floods","Dark field","Additive lag","Multiplicative lag","CNR phantom")).pack(side="left")
        self._button(demo,"Load synthetic demo",self.demo).pack(side="left",padx=5)
        rf=ttk.LabelFrame(left,text="ROIs: drag on the image, then store",padding=6);rf.pack(fill="x",pady=6)
        row=ttk.Frame(rf);row.pack(fill="x")
        for key,default in (("x",0),("y",0),("w",128),("h",128)):
            self.fields["roi_"+key]=tk.StringVar(value=str(default))
            ttk.Label(row,text=key.upper()).pack(side="left",padx=(4,2))
            ttk.Entry(row,textvariable=self.fields["roi_"+key],width=7).pack(side="left")
        self._button(row,"Set",self.set_roi).pack(side="left",padx=4)
        self._button(row,"Full",self.full_roi).pack(side="left")
        row=ttk.Frame(rf);row.pack(fill="x",pady=4)
        ttk.Combobox(row,textvariable=self.roi_slot,values=("Analysis","Signal / Inside","Background 1 / Outside","Background 2","Background 3","Background 4"),state="readonly",width=26).pack(side="left")
        self._button(row,"Store",self.store_roi).pack(side="left",padx=4)
        self._button(row,"Recall",self.recall_roi).pack(side="left")
        self._button(row,"Magic ROI",self.magic_roi).pack(side="left",padx=4)
        self.notebook=ttk.Notebook(left);self.notebook.pack(fill="both",expand=True)
        self._response_tab();self._mtf_tab();self._noise_tab();self._qc_tab();self._lag_tab();self._agd_tab()
        vpaned=ttk.Panedwindow(right,orient="vertical");vpaned.pack(fill="both",expand=True,padx=(8,0))
        imageframe=ttk.Frame(vpaned);resultframe=ttk.Frame(vpaned)
        vpaned.add(imageframe,weight=3);vpaned.add(resultframe,weight=2)
        display=ttk.Frame(imageframe);display.pack(fill="x")
        self._button(display,"Auto window",self.auto_window).pack(side="left")
        self._button(display,"Set window",self.set_window).pack(side="left",padx=4)
        ttk.Checkbutton(display,text="Use response linearization",variable=self.linearize,command=lambda:self.safe(self.redraw)).pack(side="left")
        self._button(display,"Save linearized",self.save_linearized).pack(side="left",padx=4)
        self.image_caption=tk.StringVar(value="No image loaded")
        ttk.Label(imageframe,textvariable=self.image_caption,foreground="#526173").pack(fill="x",padx=5)
        self.image_fig=Figure(figsize=(7,4),dpi=100,layout="tight");self.image_ax=self.image_fig.add_subplot(111)
        self.image_canvas=FigureCanvasTkAgg(self.image_fig,master=imageframe)
        self.image_toolbar=NavigationToolbar2Tk(self.image_canvas,imageframe,pack_toolbar=False);self.image_toolbar.pack(side="bottom",fill="x")
        self.image_canvas.get_tk_widget().pack(fill="both",expand=True)
        self.selector=RectangleSelector(self.image_ax,self.on_select,useblit=False,button=[1],minspanx=2,minspany=2,spancoords="data",interactive=True)
        self.results_notebook=ttk.Notebook(resultframe);self.results_notebook.pack(fill="both",expand=True)
        plotframe=ttk.Frame(self.results_notebook);textframe=ttk.Frame(self.results_notebook)
        self.results_notebook.add(plotframe,text="Plots");self.results_notebook.add(textframe,text="Results and details")
        resulttools=ttk.Frame(plotframe);resulttools.pack(fill="x")
        self._button(resulttools,"Export CSV / JSON",self.export).pack(side="left")
        self._button(resulttools,"Save plot",self.save_plot).pack(side="left",padx=4)
        ttk.Combobox(resulttools,textvariable=self.chart_view,values=("Curve","ESF / slit profile","LSF","2D NNPS","Signal map","SNR map"),state="readonly",width=19).pack(side="right")
        self.chart_view.trace_add("write",lambda *_:self.safe(self.plot_result))
        self.plot_fig=Figure(figsize=(7,3),dpi=100,layout="tight");self.plot_ax=self.plot_fig.add_subplot(111)
        self.plot_canvas=FigureCanvasTkAgg(self.plot_fig,master=plotframe);self.plot_canvas.get_tk_widget().pack(fill="both",expand=True)
        for canvas in (self.image_canvas,self.plot_canvas):
            canvas.get_tk_widget().bind("<Configure>",lambda _:self.root.after_idle(self._sync_canvases),add="+")
        self.root.after(150,self._sync_canvases)
        self.result_text=ScrolledText(textframe,wrap="none",font=("TkFixedFont",10));self.result_text.pack(fill="both",expand=True)
        self.result_text.insert("end","Results and acquisition settings will appear here.\nSynthetic demos are available for testing each workflow.\n")
        self._button(textframe,"Export CSV / JSON",self.export).pack(anchor="w",pady=3)
        ttk.Label(self.root,textvariable=self.status,anchor="w",padding=(12,6),wraplength=1400).grid(row=2,column=0,sticky="ew")

    def _sync_canvases(self):
        # Keep the rendering buffer inside the actual Tk widget after DPI/map changes.
        for canvas in (self.image_canvas,self.plot_canvas):
            widget=canvas.get_tk_widget();width,height=widget.winfo_width(),widget.winfo_height()
            if min(width,height)>10 and (abs(canvas.figure.bbox.width-width)>.5 or abs(canvas.figure.bbox.height-height)>.5):
                canvas.figure.set_size_inches(width/canvas.figure.dpi,height/canvas.figure.dpi,forward=False)
                canvas.draw_idle()

    def _response_tab(self):
        tab=self._tab("Response")
        self._note(tab,"Use at least four flood images at distinct measured detector air-kerma levels. Fitting uses the Analysis ROI on every loaded image, in list order.")
        self._field(tab,"Model","response_model","Linear",core.MODELS[1:])
        self._field(tab,"Air kerma list (microGy)","response_doses","",width=28)
        self._button(tab,"Fit response from all images",self.fit_response,run=True).pack(fill="x",pady=6)
        self._note(tab,"Manual calibration or identity: S = a + b*K, a + b*ln(1+c*K), a + b*ln(K), or a + b*sqrt(K). K is air kerma in microGy.")
        self._field(tab,"Calibration model","manual_model","Identity",core.MODELS)
        for k,v in (("a",0),("b",1),("c",1)):self._field(tab,k,"response_"+k,v)
        row=ttk.Frame(tab);row.pack(fill="x",pady=6)
        self._button(row,"Apply coefficients",self.manual_response).pack(side="left")
        self._button(row,"Load calibration",self.load_response).pack(side="left",padx=4)
        self._button(row,"Save calibration",self.save_response).pack(side="left")
        self.response_label=tk.StringVar(value="Current calibration: Identity")
        ttk.Label(tab,textvariable=self.response_label,wraplength=480).pack(fill="x",pady=4)

    def _mtf_tab(self):
        tab=self._tab("MTF")
        self._note(tab,"Draw the Analysis ROI around a single straight slanted edge or slit. A 1-5 degree slant with clear profiles on both sides is preferred. Physical pixel spacing is required.")
        self._field(tab,"Method","mtf_method","edge",("edge","slit"))
        self._field(tab,"Feature orientation","mtf_orientation","auto",("auto","vertical","horizontal"))
        self._field(tab,"Oversampling","mtf_os",10)
        self._field(tab,"Profile length (mm)","mtf_profile",10)
        self._field(tab,"Slit width correction (mm; 0=off)","mtf_slit",0)
        self.derivative=tk.BooleanVar(value=True)
        ttk.Checkbutton(tab,text="Correct finite-difference derivative response (edge)",variable=self.derivative).pack(anchor="w",pady=5)
        self._button(tab,"Calculate MTF",self.calculate_mtf,run=True).pack(fill="x",pady=7)
        self._button(tab,"Import MTF curve",self.import_mtf).pack(fill="x")
        self._note(tab,"Outputs: MTF, MTF50/10, edge angle, trace residual, ESF/slit profile and LSF. Frequency is limited to detector Nyquist. This implementation uses normal-distance binning and a Hamming LSF window.")

    def _noise_tab(self):
        tab=self._tab("NNPS / DQE")
        self._note(tab,"Use unprocessed, linearized flat fields acquired at the same exposure. The Analysis ROI and geometry must match across images. NNPS is area normalized and expressed in mm^2.")
        ttk.Checkbutton(tab,text="Use all loaded images for NNPS",variable=self.all_nps).pack(anchor="w")
        for label,key,default,choices in (("Fragment size (pixels)","nps_fragment",256,None),("Overlap fraction","nps_overlap",.5,None),
                                        ("Detrending","nps_detrend","mean",("mean","plane","quadratic")),("FFT window","nps_window","none",("none","hann")),
                                        ("Lines to omit near axes","nps_omit",1,None),("Lines per side to average","nps_lines",7,None)):
            self._field(tab,label,key,default,choices)
        self.reject=tk.BooleanVar(value=True)
        ttk.Checkbutton(tab,text="COQ rejection: fragment range >= absolute mean",variable=self.reject).pack(anchor="w")
        row=ttk.Frame(tab);row.pack(fill="x",pady=5)
        self._button(row,"Calculate NNPS",self.calculate_nps,run=True).pack(side="left")
        self._button(row,"Import NNPS curve",self.import_nps).pack(side="left",padx=4)
        self._field(tab,"DQE NNPS direction","nps_direction","average",("average","horizontal","vertical","radial"))
        beam=self._field(tab,"Beam (legacy COQ values)","beam","RQA 5",tuple(core.BEAMS)+("Custom",))
        beam.bind("<<ComboboxSelected>>",lambda _:self.set_beam())
        self._field(tab,"q / air kerma (1/mm^2/microGy)","q",core.BEAMS["RQA 5"])
        self._field(tab,"Detector air kerma (microGy)","dqe_kerma",2.5)
        self._field(tab,"MTF interpolation","dqe_interpolation","linear",("linear","COQ quartic"))
        self._button(tab,"Calculate DQE from stored MTF + NNPS",self.calculate_dqe,run=True).pack(fill="x",pady=5)

    def _qc_tab(self):
        tab=self._tab("QC")
        self._note(tab,"QC applies to the active image and Analysis ROI. Physical dimensions are converted separately using row and column spacing. Reported limits are legacy COQ defaults.")
        self._field(tab,"Test","qc_test","Uniformity",("Uniformity","Bad pixels","Dark analysis"))
        self._field(tab,"Uniformity ROI side (mm)","uniform_tile",30)
        self._field(tab,"Uniformity overlap","uniform_overlap",.5)
        self._field(tab,"Bad-pixel tile side (mm)","bad_tile",10)
        self._field(tab,"Bad-pixel threshold (%)","bad_threshold",20)
        self._field(tab,"Bad-pixel margin (mm)","bad_margin",15)
        self._field(tab,"Dark ROI side (mm)","dark_tile",55)
        self._field(tab,"Dark margin (mm)","dark_margin",30)
        self._button(tab,"Run selected QC test",self.calculate_qc,run=True).pack(fill="x",pady=7)
        self._note(tab,"Dark analysis reports mean and sample SD in five ROIs, matching the code's operation. Reduce ROI/margin size for small images. For CR, lag and dark analysis are unavailable; uniformity reports global values.")

    def _lag_tab(self):
        tab=self._tab("Lag / CNR")
        self._note(tab,"Store Signal / Inside and Background 1 / Outside ROIs. CNR: use 2 ROIs, or add Background 2-4 for 5 ROIs. Background pixels are pooled before computing their sample SD.")
        self._field(tab,"CNR number of ROIs","cnr_rois",2,("2","5"))
        self._button(tab,"Calculate CNR on active image",self.calculate_cnr,run=True).pack(fill="x",pady=5)
        self._button(tab,"Calculate CNR on all images",lambda:self.calculate_cnr(True),run=True).pack(fill="x")
        ttk.Separator(tab).pack(fill="x",pady=12)
        self._field(tab,"Lag mode","lag_mode","additive",("additive","multiplicative"))
        self._note(tab,"Lag uses images in list order. Additive: [test phantom, nonirradiated after]. Multiplicative: [test phantom, uniform before, uniform after]. Use exactly 2 or 3 images respectively.")
        self._button(tab,"Calculate lag",self.calculate_lag,run=True).pack(fill="x",pady=6)
        self._note(tab,"The legacy COQ lag limit is <0.5%. Preserve acquisition timing and detector correction settings when comparing measurements.")

    def _agd_tab(self):
        tab=self._tab("AGD")
        self._note(tab,"Incident air kerma must exclude backscatter. Supplied tube tables are examples. Legacy g/c polynomials are retained for comparison; use verified manual factors for your protocol.")
        self._field(tab,"PMMA thickness (mm)","agd_pmma",45)
        self._field(tab,"Anode/filter","agd_pair","Mo/Mo",("Mo/Mo","Mo/Rh","Rh/Rh","Rh/Al","W/Rh","W/Al","Custom"))
        self._field(tab,"kVp (tube lookup)","agd_kvp",28)
        self._field(tab,"HVL (mm Al)","agd_hvl",.36)
        self._field(tab,"Kerma source","agd_source","Entered incident kerma",("Entered incident kerma","Tube output"))
        self._field(tab,"Incident air kerma (mGy)","agd_kerma",5)
        self._field(tab,"mAs","agd_mas",100)
        self._field(tab,"Output (microGy/mAs)","agd_output",156)
        self._field(tab,"Output reference distance (mm)","agd_ref",500)
        self._field(tab,"Focus-to-entrance distance (mm)","agd_distance",600)
        self._field(tab,"Factor mode","agd_factors","Thickness matched COQ",("Thickness matched COQ","Original COQ index","Manual g/c/s"))
        row=ttk.Frame(tab);row.pack(fill="x",pady=3)
        for key,default in (("g",.25),("c",1.1),("s",1)):
            ttk.Label(row,text=key).pack(side="left",padx=4);self.fields["agd_"+key]=tk.StringVar(value=str(default))
            ttk.Entry(row,textvariable=self.fields["agd_"+key],width=10).pack(side="left")
        self._field(tab,"Protocol AGD limit (mGy; optional)","agd_limit","")
        row=ttk.Frame(tab);row.pack(fill="x",pady=5)
        self._button(row,"Calculate AGD",self.calculate_agd,run=True).pack(side="left")
        self._button(row,"Tube lookup",self.tube_lookup).pack(side="left",padx=3)
        self._button(row,"Load tube table",self.load_tube_table).pack(side="left")
        row=ttk.Frame(tab);row.pack(fill="x",pady=3)
        self._button(row,"DICOM acquisition fields",self.agd_metadata).pack(side="left")
        self._button(row,"Batch AGD from CSV",self.batch_agd,run=True).pack(side="left",padx=4)

    def current(self):
        if not self.images:raise ValueError("Open an image first.")
        return self.images[self.active]

    def spacing(self):
        return core.spacing2d((self.value("dy"),self.value("dx")))

    def roi(self):
        r=tuple(self.value("roi_"+k,int) for k in ("x","y","w","h"))
        core.crop(self.current().pixels,r)
        return r

    def data(self,record=None):
        record=record or self.current()
        return self.response.inverse(record.pixels) if self.linearize.get() else record.pixels

    def context(self):
        return {"modality":self.modality.get(),"images":[r.name for r in self.images],
                "image_metadata":[{"name":r.name,"shape":list(r.pixels.shape),"spacing_mm":r.spacing,"metadata":r.metadata} for r in self.images],
                "active_image":self.images[self.active].name if self.images else None,
                "spacing_mm":[self.fields[k].get() for k in ("dy","dx")],
                "analysis_roi":[self.fields["roi_"+k].get() for k in ("x","y","w","h")],
                "stored_rois":dict(self.rois),"linearization_enabled":self.linearize.get(),
                "response":asdict(self.response),"settings":{k:v.get() for k,v in self.fields.items()}}

    def invalidate(self):
        self.mtf_result=None;self.nps_result=None

    def refresh_list(self):
        self.image_list.delete(0,"end")
        for i,r in enumerate(self.images):self.image_list.insert("end",f"{i+1}. {r.name}")
        if self.images:
            self.active=min(self.active,len(self.images)-1)
            self.image_list.selection_set(self.active);self.image_list.see(self.active)
            self.select_active()
        else:
            self.image_ax.clear();self.image_ax.text(.5,.5,"Open a QC image",ha="center",va="center",transform=self.image_ax.transAxes)
            self.image_canvas.draw_idle();self.image_caption.set("No image loaded")

    def choose_image(self,_=None):
        selection=self.image_list.curselection()
        if selection:self.active=selection[0];self.safe(self.select_active)

    def select_active(self):
        rec=self.current()
        if rec.spacing:
            for k,v in zip(("dy","dx"),rec.spacing):self.fields[k].set(f"{v:g}")
        else:
            self.fields["dy"].set("");self.fields["dx"].set("")
        try:self.roi()
        except Exception:
            h,w=rec.pixels.shape
            for k,v in zip(("x","y","w","h"),((w-min(256,w))//2,(h-min(256,h))//2,min(256,w),min(256,h))):self.fields["roi_"+k].set(str(v))
        self.auto_scale=True;self.redraw()
        if not rec.spacing:self.status.set("Pixel spacing is missing. Enter measured dy and dx in mm before MTF, NNPS or QC.")

    def open_images(self):
        if self.busy:raise ValueError("Wait for the current analysis to finish before changing images.")
        paths=filedialog.askopenfilenames(title="Open original grayscale images",filetypes=[("Supported images","*.dcm *.ima *.dicom *.tif *.tiff *.npy *.png *.bmp *.jpg *.jpeg"),("All files","*")])
        if not paths:return
        errors=[]
        for p in paths:
            try:self.images.extend(core.load_images(p,self.rescale.get()))
            except Exception as e:errors.append(f"{Path(p).name}: {e}")
        self.refresh_list();self.status.set(f"{len(self.images)} images loaded. Image order and stored MTF/NNPS curves are preserved.")
        if errors:messagebox.showwarning("Some images could not be opened","\n\n".join(errors),parent=self.root)

    def clear_images(self):
        if self.busy:raise ValueError("Wait for the current calculation to finish.")
        self.images=[];self.active=0;self.rois={};self.refresh_list()

    def remove_image(self):
        self.current()
        if self.busy:raise ValueError("Wait for the current calculation to finish.")
        del self.images[self.active];self.refresh_list()

    def move_image(self,direction):
        self.current()
        if self.busy:raise ValueError("Wait for the current calculation to finish.")
        target=self.active+direction
        if 0<=target<len(self.images):
            self.images[self.active],self.images[target]=self.images[target],self.images[self.active]
            self.active=target;self.refresh_list()

    def raw_dialog(self):
        p=filedialog.askopenfilename(title="Import RAW image")
        if not p:return
        win=tk.Toplevel(self.root);win.title("RAW import");win.transient(self.root);win.grab_set()
        frame=ttk.Frame(win,padding=12);frame.pack()
        values={}
        for i,(label,key,default) in enumerate((("Rows","rows",1024),("Columns","cols",1024),("Frames","frames",1),("Offset bytes","offset",0),("Gap bytes","gap",0),("dtype","dtype","uint16"),("Byte order","byteorder","<"))):
            ttk.Label(frame,text=label).grid(row=i,column=0,sticky="w",padx=4,pady=4)
            values[key]=tk.StringVar(value=str(default))
            if key in ("dtype","byteorder"):
                opts=("uint8","uint16","int16","uint32","int32","float32","float64") if key=="dtype" else ("<",">")
                ttk.Combobox(frame,textvariable=values[key],values=opts,state="readonly").grid(row=i,column=1)
            else:ttk.Entry(frame,textvariable=values[key]).grid(row=i,column=1)
        ttk.Label(frame,text="< little endian; > big endian. Supply pixel spacing after import.").grid(row=7,columnspan=2,pady=7)
        def load():
            if self.busy:raise ValueError("Wait for the current calculation.")
            records=core.load_raw(p,**{k:v.get() for k,v in values.items()})
            self.images.extend(records);self.refresh_list();win.destroy()
        self._button(frame,"Import",load).grid(row=8,columnspan=2)

    def metadata(self):
        rec=self.current();win=tk.Toplevel(self.root);win.title("Acquisition metadata");win.geometry("700x500")
        text=ScrolledText(win,wrap="word");text.pack(fill="both",expand=True)
        text.insert("end",json.dumps({"name":rec.name,"shape":rec.pixels.shape,"spacing":rec.spacing,"metadata":rec.metadata},indent=2))
        text.configure(state="disabled")

    def on_select(self,press,release):
        if any(v is None for v in (press.xdata,press.ydata,release.xdata,release.ydata)) or not self.images:return
        h,w=self.current().pixels.shape
        x=max(0,int(np.floor(min(press.xdata,release.xdata)+.5)));y=max(0,int(np.floor(min(press.ydata,release.ydata)+.5)))
        x2=min(w,int(np.ceil(max(press.xdata,release.xdata)+.5)));y2=min(h,int(np.ceil(max(press.ydata,release.ydata)+.5)))
        if x2-x<2 or y2-y<2:return
        for k,v in zip(("x","y","w","h"),(x,y,x2-x,y2-y)):self.fields["roi_"+k].set(str(v))
        self.status.set(f"Selected ROI x={x}, y={y}, width={x2-x}, height={y2-y}. Store it to use a named signal/background ROI.")

    def set_roi(self):self.roi();self.redraw()
    def full_roi(self):
        h,w=self.current().pixels.shape
        for k,v in zip(("x","y","w","h"),(0,0,w,h)):self.fields["roi_"+k].set(str(v))
        self.redraw()
    def store_roi(self):
        self.rois[self.roi_slot.get()]=self.roi();self.redraw();self.status.set(f"Stored {self.roi_slot.get()} ROI.")
    def recall_roi(self):
        key=self.roi_slot.get()
        if key not in self.rois:raise ValueError("This ROI has not been stored.")
        for k,v in zip(("x","y","w","h"),self.rois[key]):self.fields["roi_"+k].set(str(v))
        self.redraw()
    def magic_roi(self):
        self.small_prompt("Magic ROI",(("Width in pixels","128"),("Height in pixels","128")),self._magic)
    def _magic(self,values):
        w,h=map(int,values);height,width=self.current().pixels.shape
        if w>width or h>height:raise ValueError("Magic ROI does not fit.")
        for k,v in zip(("x","y","w","h"),((width-w)//2,(height-h)//2,w,h)):self.fields["roi_"+k].set(str(v))
        self.set_roi()

    def redraw(self):
        if not self.images:return
        a=self.data();rec=self.current()
        self.selector.set_active(False);self.selector.disconnect_events();self.image_ax.clear()
        if self.auto_scale:self.display_limits=tuple(np.percentile(a,[.5,99.5]))
        lo,hi=self.display_limits
        if hi<=lo:hi=lo+1
        cmap="gray_r" if rec.metadata.get("PhotometricInterpretation")=="MONOCHROME1" else "gray"
        self.image_ax.imshow(a,cmap=cmap,origin="upper",vmin=lo,vmax=hi,interpolation="nearest")
        self.image_ax.set_xlabel("Column (pixel)");self.image_ax.set_ylabel("Row (pixel)")
        colors=["#ffcf40","#25d7f0","#ff738c","#65da91","#be8bff","#ff9955"]
        for i,(name,r) in enumerate(self.rois.items()):
            try:core.crop(a,r)
            except ValueError:continue
            x,y,w,h=r;self.image_ax.add_patch(Rectangle((x-.5,y-.5),w,h,fill=False,edgecolor=colors[i%len(colors)],linewidth=1.5))
            self.image_ax.text(x,y-1,name,fontsize=8,color=colors[i%len(colors)])
        x,y,w,h=self.roi();self.selector=RectangleSelector(self.image_ax,self.on_select,useblit=False,button=[1],minspanx=2,minspany=2,spancoords="data",interactive=True)
        self.selector.extents=(x-.5,x+w-.5,y-.5,y+h-.5)
        self.image_canvas.draw_idle()
        self.image_caption.set(f"{rec.name}  |  {a.shape[1]} x {a.shape[0]}  |  {'linearized' if self.linearize.get() else 'original signal'}  |  display [{lo:.4g}, {hi:.4g}]")

    def auto_window(self):self.auto_scale=True;self.redraw()
    def set_window(self):
        self.current();lo,hi=self.display_limits
        self.small_prompt("Display window",(("Minimum",str(lo)),("Maximum",str(hi))),self._window)
    def _window(self,values):
        lo,hi=map(float,values)
        if not np.isfinite([lo,hi]).all() or hi<=lo:raise ValueError("Window maximum must exceed minimum.")
        self.display_limits=(lo,hi);self.auto_scale=False;self.redraw()

    def small_prompt(self,title,labels,callback):
        win=tk.Toplevel(self.root);win.title(title);win.transient(self.root);win.grab_set();frame=ttk.Frame(win,padding=12);frame.pack()
        variables=[]
        for i,(label,default) in enumerate(labels):
            ttk.Label(frame,text=label).grid(row=i,column=0,sticky="w",pady=4)
            var=tk.StringVar(value=default);variables.append(var);ttk.Entry(frame,textvariable=var).grid(row=i,column=1,padx=6)
        def apply():callback([v.get() for v in variables]);win.destroy()
        self._button(frame,"Apply",apply).grid(row=len(labels),columnspan=2,pady=7)

    def save_linearized(self):
        rec=self.current();a=self.response.inverse(rec.pixels)
        p=filedialog.asksaveasfilename(title="Save linearized values",defaultextension=".tif",filetypes=[("Float TIFF","*.tif"),("NumPy array","*.npy")])
        if not p:return
        if Path(p).suffix.lower()==".npy":np.save(p,a)
        else:
            import tifffile
            tifffile.imwrite(p,a.astype(np.float32),metadata={"response":asdict(self.response),"source":rec.name,"spacing_mm":rec.spacing})
        self.status.set(f"Saved linearized data to {Path(p).name}.")

    def run(self,title,task,done=None):
        if self.busy:raise ValueError("An analysis is already running.")
        self.busy=True
        for b in self.controls:b.configure(state="disabled")
        self.status.set(f"Calculating {title}...")
        context=self.context()
        def worker():
            try:self.events.put((title,task(),None,done,context))
            except Exception as e:self.events.put((title,None,(str(e),traceback.format_exc()),None,context))
        threading.Thread(target=worker,daemon=True).start()

    def _poll(self):
        try:
            while True:
                title,result,error,done,context=self.events.get_nowait()
                self.busy=False
                for b in self.controls:b.configure(state="normal")
                if error:
                    self.status.set(error[0]);messagebox.showerror(title,error[0],parent=self.root)
                else:
                    self.safe(lambda:self._finish(title,result,context,done))
        except queue.Empty:pass
        self.root.after(100,self._poll)

    def _finish(self,title,result,context,done):
        if done:done(result)
        self.show_result(title,result,context)

    def show_result(self,kind,result,context=None):
        if kind in ("MTF","NNPS"):
            result["acquisition_context"]=context or self.context()
        self.current_result=(kind,result);self.result_context=context or self.context()
        scalars={k:v for k,v in result.items() if not isinstance(v,np.ndarray)}
        if kind=="Bad pixels":scalars["coordinates_preview_xy"]=result["coordinates_xy"][:100].tolist()
        self.result_text.delete("1.0","end")
        self.result_text.insert("end",kind+"\n"+"="*len(kind)+"\n"+json.dumps(core.serializable(scalars),indent=2)+"\n\nAcquisition and analysis context\n"+json.dumps(core.serializable(self.result_context),indent=2))
        self.chart_view.set("Curve");self.plot_result();self.status.set(f"{kind} completed. Results include settings and units; export CSV or JSON for comparison.")
        if kind=="DQE" and result.get("values_above_one"):
            self.status.set("DQE contains values >1. Check input beam q, detector kerma, linearization, MTF and NNPS acquisition conditions. Values are not clipped.")
        elif kind=="Dark analysis" or kind.startswith("AGD") or kind in ("CNR","Lag","Bad pixels"):
            self.results_notebook.select(1)

    def fit_response(self):
        if len(self.images)<4:raise ValueError("Load at least four response flood images.")
        doses=[float(v) for v in self.fields["response_doses"].get().replace(";",",").split(",") if v.strip()]
        if len(doses)!=len(self.images):raise ValueError("Enter one air-kerma value per loaded image, comma separated, in list order.")
        roi=self.roi();signals=[float(core.crop(r.pixels,roi).mean()) for r in self.images];model=self.fields["response_model"].get()
        def task():
            response=core.fit_response(doses,signals,model)
            return {"calibration":asdict(response),"rows":[{"kerma_uGy":k,"mean_signal":s,"fitted_signal":float(response.forward(k))} for k,s in zip(doses,signals)]}
        def done(result):self.response=core.Response(**result["calibration"]);self._response_changed()
        self.run("Response",task,done)

    def _response_changed(self):
        self.invalidate();self.linearize.set(self.response.model!="Identity")
        self.fields["manual_model"].set(self.response.model)
        for k in ("a","b","c"):self.fields["response_"+k].set(f"{getattr(self.response,k):.10g}")
        self.response_label.set(f"Current: {self.response.model}; a={self.response.a:.6g}, b={self.response.b:.6g}, c={self.response.c:.6g}, R²={self.response.r2}")
        if self.images:self.redraw()

    def manual_response(self):
        response=core.Response(self.fields["manual_model"].get(),*[self.value("response_"+k) for k in ("a","b","c")])
        if self.images:response.inverse(self.current().pixels)
        elif response.model!="Identity" and (abs(response.b)<1e-15 or not np.isfinite([response.a,response.b,response.c]).all()):raise ValueError("Invalid calibration parameters.")
        self.response=response;self._response_changed()
    def load_response(self):
        p=filedialog.askopenfilename(title="Load calibration JSON",filetypes=[("JSON","*.json")])
        if p:
            d=json.loads(Path(p).read_text());r=core.Response(**d.get("calibration",d));r.inverse(np.array([r.a+r.b]) if r.model=="Square root" else np.array([1.]))
            self.response=r;self._response_changed()
    def save_response(self):
        p=filedialog.asksaveasfilename(title="Save calibration",defaultextension=".json")
        if p:Path(p).write_text(json.dumps({"calibration":asdict(self.response),"units":"K=microGy, signal=input pixel units"},indent=2)+"\n")

    def calculate_mtf(self):
        a=self.data();sp=self.spacing();roi=self.roi()
        opts={"method":self.fields["mtf_method"].get(),"orientation":self.fields["mtf_orientation"].get(),"oversampling":self.value("mtf_os",int),
              "profile_mm":self.value("mtf_profile"),"derivative_correction":self.derivative.get(),"slit_width_mm":self.value("mtf_slit")}
        self.run("MTF",lambda:core.mtf(a,sp,roi,**opts),lambda r:setattr(self,"mtf_result",r))

    def calculate_nps(self):
        records=list(self.images) if self.all_nps.get() else [self.current()]
        if not records:raise ValueError("Open flat-field images first.")
        sp=self.spacing();roi=self.roi()
        if any(r.spacing and not np.allclose(r.spacing,sp) for r in records):raise ValueError("Loaded image spacings do not match entered NNPS spacing.")
        arrays=[self.data(r) for r in records]
        opts={"fragment":self.value("nps_fragment",int),"overlap":self.value("nps_overlap"),"detrend":self.fields["nps_detrend"].get(),
              "window":self.fields["nps_window"].get(),"omit":self.value("nps_omit",int),"lines":self.value("nps_lines",int),"reject_coq":self.reject.get()}
        self.run("NNPS",lambda:core.nps(arrays,sp,roi,**opts),lambda r:setattr(self,"nps_result",r))
    def import_mtf(self):
        if self.busy:raise ValueError("Wait for the current calculation.")
        p=filedialog.askopenfilename(title="Import MTF frequency/value CSV")
        if p:
            f,v=core.read_curve(p,"mtf");self.mtf_result={"frequency":f,"mtf":v,"source":"Imported "+Path(p).name};self.show_result("MTF",self.mtf_result)
    def import_nps(self):
        if self.busy:raise ValueError("Wait for the current calculation.")
        p=filedialog.askopenfilename(title="Import NNPS in mm^2 (frequency / nnps)")
        if p:
            f,v=core.read_curve(p,"average")
            curves={k:core.read_curve(p,k)[1] for k in ("average","horizontal","vertical","radial")}
            self.nps_result={"frequency":f,**curves,"source":"Imported "+Path(p).name,"units":"mm^2; confirmed by import workflow"}
            self.show_result("NNPS",self.nps_result)
    def set_beam(self):
        beam=self.fields["beam"].get()
        if beam in core.BEAMS:self.fields["q"].set(str(core.BEAMS[beam]))
    def calculate_dqe(self):
        if self.mtf_result is None or self.nps_result is None:raise ValueError("Calculate or import both MTF and NNPS first.")
        m,n=self.mtf_result,self.nps_result;direction=self.fields["nps_direction"].get()
        mc,nc=m.get("acquisition_context",{}),n.get("acquisition_context",{})
        if mc.get("linearization_enabled")!=nc.get("linearization_enabled") and "source" not in m and "source" not in n:
            raise ValueError("Stored MTF and NNPS used different linearization states. Recalculate both with the same calibration.")
        if mc.get("linearization_enabled") and mc.get("response")!=nc.get("response") and "source" not in m and "source" not in n:
            raise ValueError("Stored MTF and NNPS used different response calibrations.")
        vals=np.asarray(n[direction]);valid=np.isfinite(vals)&(vals>0)
        kerma=self.value("dqe_kerma");q=self.value("q");interpolation=self.fields["dqe_interpolation"].get()
        def task():
            result=core.dqe(m["frequency"],m["mtf"],n["frequency"][valid],vals[valid],kerma,q,interpolation)
            result["nnps_direction"]=direction;result["mtf_provenance"]={k:v for k,v in m.items() if not isinstance(v,np.ndarray)}
            result["nnps_provenance"]={k:v for k,v in n.items() if not isinstance(v,np.ndarray)}
            return result
        self.run("DQE",task)

    def calculate_qc(self):
        test=self.fields["qc_test"].get();a=self.data();sp=self.spacing();roi=self.roi()
        if test=="Dark analysis" and "CR" in self.modality.get():raise ValueError("Dark analysis is not part of the COQ CR workflow.")
        if test=="Uniformity":
            tile=self.value("uniform_tile");overlap=self.value("uniform_overlap");is_cr="CR" in self.modality.get()
            def task():
                result=core.uniformity(a,sp,roi,tile,overlap)
                if is_cr:result["LSNU_percent"]=None;result["LSNRNU_percent"]=None;result["workflow_note"]="CR: global nonuniformity only"
                return result
        elif test=="Bad pixels":
            args=[self.value(k) for k in ("bad_tile","bad_threshold","bad_margin")]
            task=lambda:core.bad_pixels(a,sp,roi,*args)
        else:
            tile=self.value("dark_tile");margin=self.value("dark_margin");task=lambda:core.dark_analysis(a,sp,roi,tile,margin)
        self.run(test,task)

    def named(self,key):
        if key not in self.rois:raise ValueError(f"Draw and store '{key}' first.")
        return self.rois[key]
    def calculate_cnr(self,batch=False):
        signal=self.named("Signal / Inside");bgs=[self.named("Background 1 / Outside")]
        if self.value("cnr_rois",int)==5:bgs += [self.named(f"Background {i}") for i in (2,3,4)]
        records=list(self.images) if batch else [self.current()]
        arrays=[self.data(r) for r in records]
        def task():
            rows=[{"image":r.name,**core.cnr(a,signal,bgs)} for a,r in zip(arrays,records)]
            return {"rows":rows,"roi_units":"zero-based pixel coordinates"}
        self.run("CNR",task)
    def calculate_lag(self):
        if "CR" in self.modality.get():raise ValueError("Lag is not part of the COQ CR workflow.")
        mode=self.fields["lag_mode"].get()
        if mode=="multiplicative" and "Angiography" in self.modality.get():raise ValueError("The COQ angiography workflow uses additive lag only.")
        inside=self.named("Signal / Inside");outside=self.named("Background 1 / Outside")
        arrays=[self.data(r) for r in self.images]
        self.run("Lag",lambda:core.lag(arrays,inside,outside,mode))

    def agd_options(self):
        mode=self.fields["agd_factors"].get();source=self.fields["agd_source"].get()
        return {"pmma_mm":self.value("agd_pmma"),"hvl_mm_al":self.value("agd_hvl"),"pair":self.fields["agd_pair"].get(),
                "incident_kerma_mGy":self.value("agd_kerma") if source=="Entered incident kerma" else None,
                "mas":self.value("agd_mas"),"output_uGy_mAs":self.value("agd_output"),
                "reference_distance_mm":self.value("agd_ref"),"entrance_distance_mm":self.value("agd_distance"),
                "factors":tuple(self.value("agd_"+k) for k in ("g","c","s")) if mode=="Manual g/c/s" else None,
                "legacy_index":mode=="Original COQ index","acceptance_limit_mGy":self.optional("agd_limit")}
    def calculate_agd(self):
        opts=self.agd_options();self.run("AGD",lambda:core.agd(**opts))
    def tube_lookup(self):
        h,o=self.tables.lookup(self.fields["agd_pair"].get(),self.value("agd_kvp"))
        self.fields["agd_hvl"].set(f"{h:g}");self.fields["agd_output"].set(f"{o:g}");self.fields["agd_ref"].set("500")
        self.status.set("Loaded HVL and output at 500 mm from the selected tube table. Verify calibration before use.")
    def load_tube_table(self):
        p=filedialog.askopenfilename(title="Load COQ two-section tube table CSV")
        if p:self.tables=core.AGDTables(p);self.tube_lookup()
    def agd_metadata(self):
        d=self.current().metadata
        mapping={"KVP":"agd_kvp","DistanceSourceToDetector":None}
        if d.get("KVP"):self.fields["agd_kvp"].set(d["KVP"])
        if d.get("ExposureInuAs"):self.fields["agd_mas"].set(str(float(d["ExposureInuAs"])/1000))
        elif d.get("Exposure"):self.fields["agd_mas"].set(d["Exposure"])
        pair="/".join(str(d.get(k,"")).upper().strip() for k in ("AnodeTargetMaterial","FilterMaterial"))
        normalize={"MOLYBDENUM":"Mo","RHODIUM":"Rh","TUNGSTEN":"W","ALUMINUM":"Al","MO":"Mo","RH":"Rh","W":"W","AL":"Al"}
        parts=pair.split("/")
        pair="/".join(normalize.get(p,p) for p in parts)
        if pair in self.tables.hvl:self.fields["agd_pair"].set(pair)
        if d.get("DistanceSourceToDetector") and d.get("BodyPartThickness"):
            self.fields["agd_distance"].set(str(float(d["DistanceSourceToDetector"])-float(d["BodyPartThickness"])))
        note="Copied kVp, mAs and available geometry. PMMA thickness and kerma require independent confirmation."
        if "EntranceDoseInmGy" in d:note+=f" DICOM EntranceDoseInmGy={d['EntranceDoseInmGy']} (derivation={d.get('EntranceDoseDerivation','unspecified')}); it has not been substituted for incident kerma."
        if "organ_dose_mGy" in d:note+=f" DICOM organ dose={d['organ_dose_mGy']} mGy (dGy converted x100)."
        messagebox.showinfo("DICOM acquisition fields",note,parent=self.root);self.status.set(note)
    def batch_agd(self):
        p=filedialog.askopenfilename(title="Batch AGD input CSV",filetypes=[("CSV","*.csv")])
        if not p:return
        with open(p,newline="",encoding="utf-8-sig") as f:rows=list(csv.DictReader(f))
        if not rows:raise ValueError("Batch CSV has no data rows.")
        def task():
            results=[]
            for i,row in enumerate(rows,2):
                try:
                    def number(k,default=None):return default if row.get(k,"").strip()=="" else float(row[k])
                    factors=tuple(number(k) for k in ("g","c","s")) if all(row.get(k,"").strip() for k in ("g","c","s")) else None
                    if any(row.get(k,"").strip() for k in ("g","c","s")) and factors is None:raise ValueError("Provide all three manual factors g/c/s.")
                    result=core.agd(number("pmma_mm"),number("hvl_mm_al"),row["anode_filter"],number("incident_kerma_mGy"),
                                    number("mas"),number("output_uGy_mAs"),number("reference_distance_mm",500),number("entrance_distance_mm"),
                                    factors=factors,legacy_index=row.get("legacy_index","").lower() in ("1","true","yes"),acceptance_limit_mGy=number("acceptance_limit_mGy"))
                    results.append({"label":row.get("label",f"row {i}"),**result})
                except Exception as e:raise ValueError(f"Batch CSV row {i}: {e}") from e
            return {"rows":results,"source":Path(p).name}
        self.run("AGD batch",task)

    def plot_result(self):
        if self.current_result is None:return
        kind,r=self.current_result;view=self.chart_view.get()
        self.plot_fig.clear();ax=self.plot_fig.add_subplot(111);self.plot_ax=ax
        if view=="2D NNPS" and "nnps2d" in r:
            arr=r["nnps2d"];fa=r["frequency_axis"];im=ax.imshow(np.log10(np.maximum(arr,1e-30)),extent=[fa.min(),fa.max(),fa.max(),fa.min()],cmap="viridis")
            self.plot_fig.colorbar(im,ax=ax,label="log10 NNPS (mm²)");ax.set_xlabel("fx (lp/mm)");ax.set_ylabel("fy (lp/mm)")
        elif view in ("Signal map","SNR map") and "signal_map" in r:
            im=ax.imshow(r["signal_map" if view=="Signal map" else "snr_map"],cmap="viridis");self.plot_fig.colorbar(im,ax=ax);ax.set_xlabel("ROI column");ax.set_ylabel("ROI row")
        elif kind=="MTF":
            if view=="ESF / slit profile" and "profile" in r:ax.plot(r["profile_x_mm"],r["profile"]);ax.set_xlabel("Normal distance (mm)");ax.set_ylabel("ESF / slit profile")
            elif view=="LSF" and "lsf" in r:ax.plot(r["lsf_x_mm"],r["lsf"]);ax.set_xlabel("Normal distance (mm)");ax.set_ylabel("LSF")
            else:
                ax.plot(r["frequency"],r["mtf"],label="MTF");ax.axhline(.5,color="#929daf",linestyle="--");ax.axhline(.1,color="#929daf",linestyle=":")
                ax.set_xlabel("Spatial frequency (lp/mm)");ax.set_ylabel("MTF");ax.set_ylim(bottom=0);ax.legend()
        elif kind=="NNPS":
            for key in ("horizontal","vertical","average","radial"):ax.semilogy(r["frequency"],r[key],label=key)
            ax.set_xlabel("Spatial frequency (lp/mm)");ax.set_ylabel("NNPS (mm²)");ax.legend(fontsize=8)
        elif kind=="DQE":ax.plot(r["frequency"],r["dqe"]);ax.set_xlabel("Spatial frequency (lp/mm)");ax.set_ylabel("DQE (fraction)");ax.set_ylim(bottom=0)
        elif kind=="Response":
            rows=r["rows"];x=[q["kerma_uGy"] for q in rows];y=[q["mean_signal"] for q in rows]
            ax.scatter(x,y,label="Measured");fine=np.linspace(min(x),max(x),200);response=core.Response(**r["calibration"]);ax.plot(fine,response.forward(fine),label="Fit")
            ax.set_xlabel("Detector air kerma (microGy)");ax.set_ylabel("Signal");ax.legend()
        elif kind=="Uniformity":
            names=["LSNU_percent","GSNU_percent","LSNRNU_percent","GSNRNU_percent"];names=[n for n in names if r[n] is not None]
            x=np.arange(len(names));ax.bar(x,[r[n] for n in names],color="#1777a8")
            ax.scatter(x,[r["limits_percent"][n] for n in names],color="#c74444",marker="_",s=160,label="Legacy COQ limit")
            ax.set_xticks(x,[n.replace("_percent","") for n in names]);ax.set_ylabel("Nonuniformity (%)");ax.legend(fontsize=8)
        elif kind=="Bad pixels":
            xy=r["coordinates_xy"]
            if len(xy):ax.scatter(xy[:,0],xy[:,1],s=2);ax.invert_yaxis()
            ax.set_xlabel("Column");ax.set_ylabel("Row");ax.set_title(f"{r['count']} flagged pixels")
        elif kind=="CNR" or kind=="AGD batch":
            key="cnr" if kind=="CNR" else "agd_mGy";rows=r["rows"]
            ax.bar(range(len(rows)),[row[key] for row in rows]);ax.set_xlabel("Image / measurement index");ax.set_ylabel("CNR" if key=="cnr" else "AGD (mGy)")
        else:ax.text(.5,.5,"Numerical results are available in Results and details.",ha="center",va="center",transform=ax.transAxes,wrap=True)
        ax.grid(alpha=.2);self.plot_canvas.draw_idle()

    def export(self):
        if not self.current_result:raise ValueError("No result to export.")
        kind,r=self.current_result
        p=filedialog.asksaveasfilename(title="Export result",initialfile="COQ_"+kind.replace(" ","_")+".csv",defaultextension=".csv",filetypes=[("CSV","*.csv"),("JSON including full arrays","*.json")])
        if p:core.export_result(p,kind,r,self.result_context);self.status.set(f"Exported {Path(p).name}.")
    def save_plot(self):
        if not self.current_result:raise ValueError("No plot to save.")
        p=filedialog.asksaveasfilename(title="Save current plot",defaultextension=".png",filetypes=[("PNG","*.png"),("PDF","*.pdf"),("SVG","*.svg")])
        if p:self.plot_fig.savefig(p,dpi=200,bbox_inches="tight")

    def demo(self):
        if self.busy:raise ValueError("Wait for the current calculation.")
        from scipy.special import erf
        self.clear_images();self.response=core.Response();self._response_changed()
        rng=np.random.default_rng(20261005);y,x=np.indices((512,512));kind=self.demo_var.get()
        records=[];self.rois={};sp=(.2,.2)
        if kind=="Slanted edge":a=1000+9000*.5*(1+erf((x-256-.055*(y-256))/(np.sqrt(2)*1.1)));arrays=[a+rng.normal(0,2,a.shape)]
        elif kind=="Slanted slit":a=1000+9000*np.exp(-.5*((x-256-.055*(y-256))/1.1)**2);arrays=[a+rng.normal(0,2,a.shape)]
        elif kind=="Flat fields":arrays=[10000+rng.normal(0,100,(512,512)) for _ in range(4)]
        elif kind=="Response floods":
            doses=[.5,1,2,4,8];arrays=[500+1500*k+rng.normal(0,10,(512,512)) for k in doses]
            self.fields["response_doses"].set(",".join(map(str,doses)))
        elif kind=="Dark field":arrays=[500+rng.normal(0,8,(512,512))]
        elif kind=="CNR phantom":
            a=10000+rng.normal(0,100,(512,512));a[200:250,200:250]-=500;arrays=[a]
            self.rois={"Signal / Inside":(205,205,40,40),"Background 1 / Outside":(280,205,40,40)}
        else:
            test=np.where(x<256,1000.,10000.);after=np.where(x<256,502.,500.)
            if kind=="Additive lag":arrays=[test,after];self.fields["lag_mode"].set("additive")
            else:arrays=[test,np.full(x.shape,10000.),np.where(x<256,9970.,10000.)];self.fields["lag_mode"].set("multiplicative")
            self.rois={"Signal / Inside":(100,200,40,40),"Background 1 / Outside":(350,200,40,40)}
        self.images=[core.ImageRecord(f"SYNTHETIC {kind} {i+1}",a,sp,{"synthetic":True}) for i,a in enumerate(arrays)]
        self.refresh_list();self.full_roi()
        if "Slanted" in kind:
            for k,v in zip(("x","y","w","h"),(128,128,256,256)):self.fields["roi_"+k].set(str(v))
            self.fields["mtf_method"].set("slit" if "slit" in kind else "edge");self.notebook.select(1)
        elif kind=="Flat fields":self.notebook.select(2)
        elif kind=="Response floods":self.notebook.select(0)
        elif kind=="Dark field":self.fields["dark_tile"].set("10");self.fields["dark_margin"].set("5");self.fields["qc_test"].set("Dark analysis");self.notebook.select(3)
        else:self.notebook.select(4)
        self.redraw();self.status.set(f"Loaded {kind} demonstration. Synthetic data are for software testing.")

    def help(self):
        win=tk.Toplevel(self.root);win.title("COQ Python help");win.geometry("850x680")
        text=ScrolledText(win,wrap="word",font=("TkDefaultFont",11));text.pack(fill="both",expand=True)
        p=Path(__file__).with_name("README.md")
        text.insert("end",p.read_text(encoding="utf-8") if p.exists() else "Read README.md for workflows, equations and differences from the original plugin.")
        text.configure(state="disabled")

    def close(self):
        self.root.destroy()


def main():
    parser=argparse.ArgumentParser(description="Standalone COQ detector QA GUI")
    parser.add_argument("images",nargs="*",help="Image files to open at startup")
    args=parser.parse_args()
    root=tk.Tk();app=COQApp(root)
    if args.images:
        for p in args.images:
            try:app.images.extend(core.load_images(p))
            except Exception as e:messagebox.showerror("Open image",str(e),parent=root)
        app.refresh_list()
    root.mainloop()


if __name__=="__main__":main()
