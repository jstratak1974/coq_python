"""Exercise GUI callbacks on a desktop/display. Run python tests/gui_smoke_test.py.

Dialogs are replaced with temporary paths during this test. No patient data used.
"""
import json
import os
import sys
import tempfile
import time
from pathlib import Path
import tkinter as tk

sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from coq_python import COQApp, filedialog, messagebox
import coq_core as c


def exercise(screenshot=None):
    root=tk.Tk();app=COQApp(root);errors=[];checks=[]
    messagebox.showerror=lambda title,message,**kw:errors.append(f"{title}: {message}")
    messagebox.showwarning=lambda *a,**kw:None
    def wait():
        deadline=time.monotonic()+20
        while app.busy and time.monotonic()<deadline:
            root.update();time.sleep(.01)
        root.update()
        assert not app.busy,"GUI calculation timed out"
        assert not errors,"; ".join(errors)
    def demo(name):app.demo_var.set(name);app.demo();root.update()
    try:
        demo("Slanted edge");app.calculate_mtf();wait();checks.append("edge MTF callback")
        stored_mtf=app.mtf_result
        # Loading/clearing images must preserve MTF for a later flat-field acquisition.
        app.clear_images();assert app.mtf_result is stored_mtf
        demo("Slanted slit");app.calculate_mtf();wait();checks.append("slit MTF callback")
        stored_mtf=app.mtf_result
        demo("Flat fields");app.mtf_result=stored_mtf
        app.calculate_nps();wait();checks.append("NNPS callback")
        app.calculate_dqe();wait();checks.append("DQE callback")
        with tempfile.TemporaryDirectory() as folder:
            p=Path(folder)
            for kind,r in (("MTF",app.mtf_result),("NNPS",app.nps_result)):
                c.export_result(p/(kind+".csv"),kind,r)
            filedialog.askopenfilename=lambda **kw:str(p/"MTF.csv")
            app.import_mtf()
            filedialog.askopenfilename=lambda **kw:str(p/"NNPS.csv")
            app.import_nps();app.calculate_dqe();wait();checks.append("MTF/NNPS CSV import and DQE")
            filedialog.asksaveasfilename=lambda **kw:str(p/"result.json")
            app.export();assert json.loads((p/"result.json").read_text())["analysis"]=="DQE"
            checks.append("JSON export")
            app.show_result("NNPS",app.nps_result)
            filedialog.asksaveasfilename=lambda **kw:str(p/"result.csv")
            app.export();assert (p/"result.csv").stat().st_size>100
            app.chart_view.set("2D NNPS");root.update()
            filedialog.asksaveasfilename=lambda **kw:str(p/"plot.png")
            app.save_plot();assert (p/"plot.png").stat().st_size>100
            checks.append("CSV and plot export")
        for test in ("Uniformity","Bad pixels"):
            app.fields["qc_test"].set(test);app.calculate_qc();wait();checks.append(test+" callback")
        demo("Dark field");app.calculate_qc();wait();checks.append("dark callback")
        demo("Additive lag");app.calculate_lag();wait();checks.append("additive lag callback")
        demo("Multiplicative lag");app.calculate_lag();wait();checks.append("multiplicative lag callback")
        demo("CNR phantom");app.calculate_cnr();wait();checks.append("CNR callback")
        app.calculate_cnr(True);wait();checks.append("batch CNR callback")
        demo("Response floods");app.fit_response();wait();assert app.response.r2>.9999
        checks.append("response fitting and automatic linearization")
        app.calculate_agd();wait();checks.append("AGD callback")
        app.fields["agd_source"].set("Tube output");app.tube_lookup();app.calculate_agd();wait();checks.append("tube-output AGD callback")
        filedialog.askopenfilename=lambda **kw:str(Path(__file__).resolve().parents[1]/"agd_batch_template.csv")
        app.batch_agd();wait();assert len(app.current_result[1]["rows"])==3;checks.append("AGD batch CSV callback")
        demo("Slanted edge");app.calculate_mtf();wait()
        app.results_notebook.select(0);app.root.geometry("1480x940+0+0");root.update()
        if screenshot:
            from PIL import ImageGrab
            ImageGrab.grab(xdisplay=os.environ.get("DISPLAY")).save(screenshot)
        return {"checks":checks,"passed":len(checks),"errors":errors}
    finally:root.destroy()


if __name__=="__main__":print(json.dumps(exercise(),indent=2))
