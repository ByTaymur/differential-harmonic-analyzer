"""
 Profesyonel Harmonik Analizör - Dual Channel Analyzer
 ======================================================
 - Tek CSV dosyasından CH1 ve CH2 okuma
 - PNG görüntüden dalga formu çıkarma
 - Iki kanalda bağımsız akım ölçümü
 - Her kanal için ayrı ratio (A/V) dönüşümü
 - Akımları ayrı veya üst üste görme
 - IEC 61000-3-2 harmonik limitleri
 - Power Factor, TDD, Phase analysis
 - Batch processing
 - Export to CSV/Excel

 Akım Probu Dönüşümü: 5A -> 0.25V demek ratio = 20 A/V
"""

import tkinter as tk
from tkinter import ttk, filedialog, messagebox, simpledialog
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.widgets import Cursor
import numpy as np
from scipy.fft import fft, fftfreq
from scipy.signal import butter, filtfilt, savgol_filter, find_peaks
from scipy.ndimage import gaussian_filter1d
from datetime import datetime
import os
import json
from collections import defaultdict
import threading

# IEC 61000-3-2 CLASS A LIMITLERI (Amper)
IEC_CLASS_A_LIMITS = {
    2: 1.0800, 3: 2.3000, 4: 0.4300, 5: 1.1400, 6: 0.3000,
    7: 0.7700, 8: 0.2300, 9: 0.4000, 10: 0.1840, 11: 0.3300,
    12: 0.1533, 13: 0.2100, 14: 0.1314, 15: 0.1500, 16: 0.1150,
    17: 0.1324, 18: 0.1022, 19: 0.1184, 20: 0.0920, 21: 0.1071,
    22: 0.0836, 23: 0.0978, 24: 0.0767, 25: 0.0900, 26: 0.0708,
    27: 0.0833, 28: 0.0657, 29: 0.0776, 30: 0.0613, 31: 0.0726,
    32: 0.0575, 33: 0.0682, 34: 0.0541, 35: 0.0643, 36: 0.0511,
    37: 0.0608, 38: 0.0484, 39: 0.0577, 40: 0.0460,
}

# Preset dönüşüm oranları (A/V)
RATIO_PRESETS = {
    "5A->0.25V (20 A/V)": 20.0,
    "10A->1V (10 A/V)": 10.0,
    "1A->0.1V (10 A/V)": 10.0,
    "1A->1V (1 A/V)": 1.0,
    "100mV/A (10 A/V)": 10.0,
    "50mV/A (20 A/V)": 20.0,
    "Manual": None
}

# Analiz presetleri
ANALYSIS_PRESETS = {
    "IEC61000-3-2 Class A": {
        "harmonics": 40,
        "fundamental_range": (45, 65),
        "thd_limit": 100,
        "limits": IEC_CLASS_A_LIMITS
    },
    "Hızlı Analiz": {
        "harmonics": 20,
        "fundamental_range": (45, 65),
        "thd_limit": 100,
        "limits": {k: IEC_CLASS_A_LIMITS[k] for k in range(2, 21)}
    },
    "Geniş Bant": {
        "harmonics": 50,
        "fundamental_range": (45, 65),
        "thd_limit": 100,
        "limits": IEC_CLASS_A_LIMITS
    }
}


class ImageWaveformExtractor:
    """PNG görüntüden dalga formu çıkarma sınıfı"""
    
    def __init__(self):
        self.calibration = {
            'x0': 50, 'x1': 750,  # Grid sınırları
            'y0': 50, 'y1': 550,
            'time_scale': 0.02,  # 20ms tam skala
            'volt_scale': 1.0,   # 1V tam skala
            'time_unit': 's',
            'volt_unit': 'V'
        }
    
    def extract_waveform(self, image_path, calibration=None):
        """Görüntüden dalga formu çıkar"""
        try:
            from PIL import Image, ImageOps
            import cv2
            import numpy as np
            
            if calibration:
                self.calibration.update(calibration)
            
            # Görüntüyü aç
            img = Image.open(image_path)
            img_gray = ImageOps.grayscale(img)
            img_array = np.array(img_gray)
            
            # Threshold ile waveform detection
            _, binary = cv2.threshold(img_array, 127, 255, cv2.THRESH_BINARY)
            
            # Morfolojik işlemler
            kernel = np.ones((3, 3), np.uint8)
            binary = cv2.morphologyEx(binary, cv2.MORPH_CLOSE, kernel)
            
            # Grid alanını belirle
            x_start, x_end = self.calibration['x0'], self.calibration['x1']
            y_start, y_end = self.calibration['y0'], self.calibration['y1']
            
            # Her x için y değerini bul
            waveform_y = []
            waveform_x = []
            
            grid_width = x_end - x_start
            grid_height = y_end - y_start
            
            for i in range(grid_width):
                col = x_start + i
                col_data = binary[y_start:y_end, col]
                
                # En üst non-zero pixel'i bul (waveform tipik olarak grid üstünde)
                nonzero_indices = np.where(col_data < 255)[0]
                
                if len(nonzero_indices) > 0:
                    # Tipik olarak grid ortasından yukarı çizim
                    # En yakın sinyal pixel'ini bul
                    y_idx = y_start + (y_end - y_start) - nonzero_indices[0]
                    waveform_y.append(y_idx)
                    waveform_x.append(i)
            
            if len(waveform_x) < 10:
                return None
            
            # Normalize et
            time = np.array(waveform_x) / len(waveform_x) * self.calibration['time_scale']
            signal = (self.calibration['y1'] - np.array(waveform_y)) / (self.calibration['y1'] - self.calibration['y0'])
            signal = (signal - 0.5) * 2 * self.calibration['volt_scale']
            
            # Interpolasyon ile düzleştir
            from scipy.interpolate import interp1d
            f = interp1d(time, signal, kind='cubic', fill_value='extrapolate')
            time_smooth = np.linspace(time.min(), time.max(), len(time) * 2)
            signal_smooth = f(time_smooth)
            
            sample_rate = len(time) / (time.max() - time.min()) if time.max() > time.min() else 10000
            
            return {
                'time': time_smooth,
                'signal': signal_smooth,
                'sample_rate': sample_rate,
                'calibration': self.calibration.copy()
            }
            
        except ImportError:
            messagebox.showwarning("Eksik Kütüphane", 
                "PIL ve OpenCV kurulu değil. pip install pillow opencv-python")
            return None
        except Exception as e:
            messagebox.showerror("Hata", f"Görüntü işleme hatası: {str(e)}")
            return None


class HarmonicAnalyzer:
    """Profesyonel Harmonik Analiz Sınıfı - Labaratuvar Cihazı Uyumlu"""
    
    def __init__(self):
        self.iec_limits = IEC_CLASS_A_LIMITS
    
    def calculate_all_metrics(self, signal, sample_rate, fundamental_freq=None, num_harmonics=40):
        """Tüm metrikleri hesapla - harmonik_simple.py ve iec_harmonic_analyzer.py ile uyumlu"""
        # DC offset kaldır (opsiyonel - harmonik analiz için önemli değil ama temiz veri için)
        signal = signal - np.mean(signal)
        
        # Temel frekans bul - harmonik_simple.py ile aynı yöntem
        if fundamental_freq is None:
            fundamental_freq = self.find_fundamental(signal, sample_rate)
        
        # Harmonik analizi - iec_harmonic_analyzer.py yöntemi ile aynı
        harmonics = self.calculate_harmonics_standard(signal, sample_rate, fundamental_freq, num_harmonics)
        
        # THD hesapla
        thd = self.calculate_thd(harmonics)
        
        # TDD hesapla
        tdd = self.calculate_tdd(harmonics)
        
        # RMS hesapla
        rms = np.sqrt(np.mean(signal**2))
        
        # Peak değerleri
        ipk = np.max(np.abs(signal))
        
        # Crest Factor
        cf = ipk / rms if rms > 0 else 0
        
        # Power Factor
        pf = self.calculate_power_factor(signal, sample_rate, fundamental_freq)
        
        return {
            'fundamental': fundamental_freq,
            'harmonics': harmonics,
            'thd': thd,
            'tdd': tdd,
            'rms': rms,
            'ipk': ipk,
            'cf': cf,
            'pf': pf,
            'ff': cf,
            'passed': self.check_iec_compliance(harmonics),
            'failed': [h for h in harmonics if h['status'] == 'FAIL']
        }
    
    def find_fundamental(self, signal, sample_rate):
        """Temel frekansı bul - harmonik_simple.py ile aynı"""
        n = len(signal)
        yf = np.abs(fft(signal))[:n//2]
        xf = fftfreq(n, 1/sample_rate)[:n//2]
        
        # 45-65 Hz arası ara
        mask = (xf >= 45) & (xf <= 65)
        idx = np.where(mask)[0]
        
        if len(idx) > 0:
            peak_idx = idx[np.argmax(yf[idx])]
            return xf[peak_idx]
        return 50.0
    
    def calculate_harmonics_standard(self, signal, sample_rate, fundamental, num_harmonics=40):
        """Standart harmonik hesaplama - iec_harmonic_analyzer.py ile aynı"""
        n = len(signal)
        
        # Tam FFT - pencereleme YOK (lab cihazları gibi)
        yf_full = fft(signal)
        xf_full = fftfreq(n, 1/sample_rate)
        
        positive_mask = xf_full >= 0
        xf_pos = xf_full[positive_mask]
        
        # Genlik hesaplama: 2/n ölçekleme (tepe genlik için)
        yf_pos = np.abs(yf_full[positive_mask]) * 2 / n
        
        harmonics = []
        for h in range(1, num_harmonics + 1):
            target_freq = h * fundamental
            
            # Hedef frekansın indeksini bul
            idx = np.argmin(np.abs(xf_pos - target_freq))
            
            # ±3 bin lokal arama (harmonik_simple.py ile aynı)
            search_start = max(0, idx - 3)
            search_end = min(len(yf_pos), idx + 4)
            local_max_idx = search_start + np.argmax(yf_pos[search_start:search_end])
            
            amplitude = yf_pos[local_max_idx]
            
            # Faz hesabı
            phase = np.angle(yf_full[positive_mask][local_max_idx]) * 180 / np.pi
            
            # Limit kontrolü
            limit = self.iec_limits.get(h, 0) if h > 1 else 0
            percent = (amplitude / limit * 100) if limit > 0 else 0
            
            harmonics.append({
                'harmonic': h,
                'frequency': target_freq,
                'amplitude': amplitude,
                'phase': phase,
                'limit': limit,
                'percent': percent,
                'status': 'FUND' if h == 1 else ('FAIL' if percent > 100 else 'PASS')
            })
        
        return harmonics
    
    def calculate_thd(self, harmonics):
        """THD hesapla - harmonik_simple.py ile aynı"""
        fundamental = harmonics[0]['amplitude']
        if fundamental == 0:
            return 0
        sum_squares = sum(h['amplitude']**2 for h in harmonics[1:41])
        return np.sqrt(sum_squares) / fundamental * 100
    
    def calculate_tdd(self, harmonics, fundamental_rms=None):
        """TDD hesapla"""
        if fundamental_rms is None:
            fundamental_rms = harmonics[0]['amplitude']
        if fundamental_rms == 0:
            return 0
        sum_squares = sum(h['amplitude']**2 for h in harmonics[1:41])
        return np.sqrt(sum_squares) / fundamental_rms * 100
    
    def calculate_power_factor(self, signal, sample_rate, fundamental):
        """Güç faktörü hesapla - iec_harmonic_analyzer.py yöntemi"""
        n = len(signal)
        
        # Basit yaklaşım: temel harmonik genliği / RMS
        yf_full = fft(signal)
        xf_full = fftfreq(n, 1/sample_rate)
        
        idx = np.argmin(np.abs(xf_full - fundamental))
        fundamental_amplitude = np.abs(yf_full[idx]) * 2 / n
        
        rms = np.sqrt(np.mean(signal**2))
        
        if rms > 0 and fundamental_amplitude > 0:
            # Displacement factor yaklaşık
            displacement = fundamental_amplitude / rms
            return min(1.0, displacement)
        
        return 1.0
    
    def check_iec_compliance(self, harmonics):
        """IEC uyumluluğunu kontrol et"""
        for h in harmonics[1:]:
            if h['status'] == 'FAIL':
                return False
        return True
        sum_squares = sum(h['amplitude']**2 for h in harmonics[1:])
        return np.sqrt(sum_squares) / fundamental_rms * 100
    
    def calculate_power_factor(self, signal, sample_rate, fundamental):
        """Güç faktörü hesapla"""
        n = len(signal)
        yf = fft(signal)
        xf = fftfreq(n, 1/sample_rate)
        
        # Temel frekans indeksini bul
        idx = np.argmin(np.abs(xf - fundamental))
        if idx < len(yf):
            fundamental_amplitude = np.abs(yf[idx]) * 2 / n
            # Basit PF hesabı
            return min(1.0, fundamental_amplitude / (np.sqrt(np.mean(signal**2)) + 0.0001))
        return 1.0
    
    def check_iec_compliance(self, harmonics):
        """IEC uyumluluğunu kontrol et"""
        for h in harmonics[1:]:  # Temel hariç
            if h['status'] == 'FAIL':
                return False
        return True


class DualCurrentAnalyzer:
    def __init__(self, root):
        self.root = root
        self.root.title("🔬 Profesyonel Harmonik Analizör - IEC 61000-3-2")
        self.root.geometry("1800x1100")
        self.root.configure(bg='#1a1a2e')
        
        self.data = {}
        self.results = {}
        self.analyzer = HarmonicAnalyzer()
        self.image_extractor = ImageWaveformExtractor()
        
        # Cursor state
        self.cursor1_pos = None
        self.cursor2_pos = None
        self.active_cursor = None
        
        # Batch processing
        self.batch_files = []
        self.batch_index = 0
        
        self.setup_styles()
        self.setup_ui()
        self.setup_shortcuts()
    
    def setup_styles(self):
        """GUI stilleri"""
        self.style = ttk.Style()
        self.style.theme_use('clam')
        
        bg_dark = '#1a1a2e'
        bg_card = '#16213e'
        bg_input = '#0f3460'
        
        self.style.configure('Dark.TFrame', background=bg_dark)
        self.style.configure('Card.TFrame', background=bg_card)
        self.style.configure('Dark.TLabel', background=bg_dark, foreground='#e8e8e8', font=('Segoe UI', 10))
        self.style.configure('Title.TLabel', background=bg_dark, foreground='#00d4ff', font=('Segoe UI', 14, 'bold'))
        self.style.configure('CH1.TLabel', background=bg_dark, foreground='#00d4ff', font=('Segoe UI', 11, 'bold'))
        self.style.configure('CH2.TLabel', background=bg_dark, foreground='#ff8844', font=('Segoe UI', 11, 'bold'))
        self.style.configure('Pass.TLabel', background=bg_dark, foreground='#00ff88', font=('Segoe UI', 14, 'bold'))
        self.style.configure('Fail.TLabel', background=bg_dark, foreground='#ff4444', font=('Segoe UI', 14, 'bold'))
        self.style.configure('Status.TLabel', background=bg_dark, foreground='#888888', font=('Segoe UI', 9))
        self.style.configure('TLabelframe', background=bg_dark, foreground='white')
        self.style.configure('TLabelframe.Label', background=bg_dark, foreground='white')
    
    def setup_ui(self):
        """Ana arayüz"""
        main_frame = ttk.Frame(self.root, style='Dark.TFrame', padding="5")
        main_frame.pack(fill=tk.BOTH, expand=True)
        
        # === ÜST PANEL: Sekmeler ===
        notebook = ttk.Notebook(main_frame)
        notebook.pack(fill=tk.X, pady=(0, 5))
        
        # Ana sekme
        self.main_tab = ttk.Frame(notebook, style='Dark.TFrame')
        notebook.add(self.main_tab, text="📊 Ana Analiz")
        
        # Batch sekmesi
        self.batch_tab = ttk.Frame(notebook, style='Dark.TFrame')
        notebook.add(self.batch_tab, text="📁 Batch İşlem")
        
        # Rapor sekmesi
        self.report_tab = ttk.Frame(notebook, style='Dark.TFrame')
        notebook.add(self.report_tab, text="📋 Rapor")
        
        self.setup_main_tab()
        self.setup_batch_tab()
        self.setup_report_tab()
        
        # === DURUM ÇUBUĞU ===
        self.status_bar = ttk.Label(main_frame, text="Hazır | F1: Yardım | Ctrl+O: Dosya | Ctrl+I: Görüntü",
                                    style='Status.TLabel', anchor='w')
        self.status_bar.pack(fill=tk.X, pady=(5, 0))
    
    def setup_main_tab(self):
        """Ana sekme arayüzü"""
        # Sol panel - Kontroller
        control_panel = ttk.Frame(self.main_tab, style='Card.TFrame', padding="10")
        control_panel.pack(side=tk.LEFT, fill=tk.Y, padx=(0, 5))
        
        # Dosya seçimi
        self.create_file_section(control_panel)
        
        # Kanal ayarları
        self.create_channel_section(control_panel, 'CH1', '#00d4ff')
        self.create_channel_section(control_panel, 'CH2', '#ff8844')
        
        # Analiz ayarları
        self.create_analysis_section(control_panel)
        
        # İşlem butonları
        self.create_action_section(control_panel)
        
        # Sağ panel - Grafikler
        graph_panel = ttk.Frame(self.main_tab, style='Dark.TFrame')
        graph_panel.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        
        # Sonuç banner
        self.result_banner = ttk.Label(graph_panel, text="CSV veya PNG dosyası yükleyin",
                                       style='Dark.TLabel', font=('Consolas', 12))
        self.result_banner.pack(fill=tk.X, pady=(0, 5))
        
        # Grafik alanı
        self.fig = plt.figure(figsize=(12, 8))
        self.fig.patch.set_facecolor('#1a1a2e')
        
        self.canvas = FigureCanvasTkAgg(self.fig, master=graph_panel)
        self.canvas.draw()
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        
        # Özel Araç Çubuğu
        self.create_custom_toolbar(graph_panel)
        
        # Cursor event
        self.canvas.mpl_connect('button_press_event', self.on_canvas_click)
        self.canvas.mpl_connect('motion_notify_event', self.on_mouse_move)
    
    def create_custom_toolbar(self, parent):
        """Özel araç çubuğu - Zoom, Pan, Referans"""
        toolbar_frame = ttk.Frame(parent)
        toolbar_frame.pack(fill=tk.X, pady=(5, 0))
        
        # Zoom butonları
        ttk.Button(toolbar_frame, text="🔍+", width=5, command=self.zoom_in).pack(side=tk.LEFT, padx=2)
        ttk.Button(toolbar_frame, text="🔍-", width=5, command=self.zoom_out).pack(side=tk.LEFT, padx=2)
        ttk.Button(toolbar_frame, text="🏠", width=5, command=self.zoom_reset).pack(side=tk.LEFT, padx=2)
        
        ttk.Separator(toolbar_frame, orient='vertical').pack(side=tk.LEFT, fill=tk.Y, padx=10)
        
        # Pan butonları
        ttk.Button(toolbar_frame, text="✋", width=5, command=self.toggle_pan).pack(side=tk.LEFT, padx=2)
        ttk.Button(toolbar_frame, text="↔️", width=5, command=self.pan_left).pack(side=tk.LEFT, padx=2)
        ttk.Button(toolbar_frame, text="↔️", width=5, command=self.pan_right).pack(side=tk.LEFT, padx=2)
        
        ttk.Separator(toolbar_frame, orient='vertical').pack(side=tk.LEFT, fill=tk.Y, padx=10)
        
        # Referans çizgisi
        ttk.Button(toolbar_frame, text="📏 Set Ref", width=10, command=self.set_reference_line).pack(side=tk.LEFT, padx=2)
        ttk.Button(toolbar_frame, text="❌ Clear Ref", width=10, command=self.clear_reference_line).pack(side=tk.LEFT, padx=2)
        
        # Zoom bölgesi seç
        ttk.Button(toolbar_frame, text="📐 Zoom Box", width=12, command=self.enable_zoom_box).pack(side=tk.LEFT, padx=2)
        
        ttk.Separator(toolbar_frame, orient='vertical').pack(side=tk.LEFT, fill=tk.Y, padx=10)
        
        # Coordinates display
        self.coord_label = ttk.Label(toolbar_frame, text="x: --, y: --", style='Status.TLabel', width=25)
        self.coord_label.pack(side=tk.LEFT, padx=5)
        
        # Reference info
        self.ref_label = ttk.Label(toolbar_frame, text="", style='Status.TLabel', foreground='#00ff88')
        self.ref_label.pack(side=tk.LEFT, padx=5)
    
    def zoom_in(self):
        """Zoom in - merkezi yaklaştır"""
        for ax in self.fig.axes:
            xlim = ax.get_xlim()
            ylim = ax.get_ylim()
            x_center = (xlim[0] + xlim[1]) / 2
            y_center = (ylim[0] + ylim[1]) / 2
            x_range = (xlim[1] - xlim[0]) * 0.4
            y_range = (ylim[1] - ylim[0]) * 0.4
            ax.set_xlim(x_center - x_range, x_center + x_range)
            ax.set_ylim(y_center - y_range, y_center + y_range)
        self.canvas.draw()
    
    def zoom_out(self):
        """Zoom out - merkezi uzaklaştır"""
        for ax in self.fig.axes:
            xlim = ax.get_xlim()
            ylim = ax.get_ylim()
            x_center = (xlim[0] + xlim[1]) / 2
            y_center = (ylim[0] + ylim[1]) / 2
            x_range = (xlim[1] - xlim[0]) * 1.5
            y_range = (ylim[1] - ylim[0]) * 1.5
            ax.set_xlim(x_center - x_range, x_center + x_range)
            ax.set_ylim(y_center - y_range, y_center + y_range)
        self.canvas.draw()
    
    def zoom_reset(self):
        """Tüm zoom/pan'ı sıfırla"""
        self.fig.tight_layout()
        self.canvas.draw()
        self.update_plots()
    
    def toggle_pan(self):
        """Pan modunu aç/kapa"""
        # Pan modu için sağ tık ile sürükleme
        self.canvas.mpl_connect('button_press_event', self.on_pan_click)
        self.status_bar.config(text="Pan: Sürüklemek için farenin sağ tuşunu basılı tutun")
    
    def pan_left(self):
        """Sola kaydır"""
        for ax in self.fig.axes:
            xlim = ax.get_xlim()
            shift = (xlim[1] - xlim[0]) * 0.2
            ax.set_xlim(xlim[0] - shift, xlim[1] - shift)
        self.canvas.draw()
    
    def pan_right(self):
        """Sağa kaydır"""
        for ax in self.fig.axes:
            xlim = ax.get_xlim()
            shift = (xlim[1] - xlim[0]) * 0.2
            ax.set_xlim(xlim[0] + shift, xlim[1] + shift)
        self.canvas.draw()
    
    def on_pan_click(self, event):
        """Pan için fare olayı"""
        if event.button == 3 and event.inaxes:  # Sağ tık
            self._pan_start = (event.xdata, event.ydata)
            self._pan_axes = event.inaxes
            self.canvas.mpl_connect('motion_notify_event', self.on_pan_drag)
            self.canvas.mpl_connect('button_release_event', self.on_pan_release)
    
    def on_pan_drag(self, event):
        """Pan sürükleme"""
        if hasattr(self, '_pan_start') and event.inaxes == self._pan_axes:
            dx = self._pan_start[0] - event.xdata
            dy = self._pan_start[1] - event.ydata
            ax = self._pan_axes
            xlim = ax.get_xlim()
            ylim = ax.get_ylim()
            ax.set_xlim(xlim[0] + dx, xlim[1] + dx)
            ax.set_ylim(ylim[0] + dy, ylim[1] + dy)
            self.canvas.draw_idle()
    
    def on_pan_release(self, event):
        """Pan bitti"""
        self.canvas.mpl_disconnect(self._pan_drag_id)
        self.canvas.mpl_disconnect(self._pan_release_id)
    
    def enable_zoom_box(self):
        """Kutu ile zoom seçimi"""
        self.canvas.mpl_connect('button_press_event', self.on_zoom_box_start)
        self.status_bar.config(text="Zoom Box: İlk köşeyi seçin (sol tık)")
    
    def on_zoom_box_start(self, event):
        """Zoom box başlangıç"""
        if event.button == 1 and event.inaxes:  # Sol tık
            self._zoom_start = (event.xdata, event.ydata)
            self._zoom_axes = event.inaxes
            self.canvas.mpl_connect('motion_notify_event', self.on_zoom_box_drag)
            self.canvas.mpl_connect('button_release_event', self.on_zoom_box_end)
    
    def on_zoom_box_drag(self, event):
        """Zoom box sürükleme - dikdörtgen çiz"""
        if hasattr(self, '_zoom_start') and event.inaxes == self._zoom_axes:
            ax = self._zoom_axes
            # Mevcut limits
            xlim = ax.get_xlim()
            ylim = ax.get_ylim()
            
            # Geçici dikdörtgen çiz (silme)
            self.fig.canvas.draw()
    
    def on_zoom_box_end(self, event):
        """Zoom box bitiş"""
        if hasattr(self, '_zoom_start') and event.inaxes == self._zoom_axes:
            x0, y0 = self._zoom_start
            x1, y1 = event.xdata, event.ydata
            
            ax = self._zoom_axes
            xlim = ax.get_xlim()
            ylim = ax.get_ylim()
            
            # Yeni sınırları ayarla
            new_xlim = (min(x0, x1), max(x0, x1))
            new_ylim = (min(y0, y1), max(y0, y1))
            
            if new_xlim[0] != new_xlim[1] and new_ylim[0] != new_ylim[1]:
                ax.set_xlim(new_xlim)
                ax.set_ylim(new_ylim)
                self.canvas.draw()
            
            self.canvas.mpl_disconnect(self._zoom_drag_id)
            self.canvas.mpl_disconnect(self._zoom_end_id)
            self.status_bar.config(text="Zoom uygulandı")
    
    def set_reference_line(self):
        """Referans çizgisi ayarla (tıklanan nokta)"""
        if self.cursor1_pos:
            self.ref_x = self.cursor1_pos[0]
            self.ref_y = self.cursor1_pos[1]
            self.ref_label.config(text=f"REF: x={self.ref_x:.4f}, y={self.ref_y:.4f}")
            
            # Çizgi ekle
            for ax in self.fig.axes:
                ax.axvline(self.ref_x, color='yellow', linestyle='--', linewidth=1, alpha=0.7)
                ax.axhline(self.ref_y, color='yellow', linestyle='--', linewidth=1, alpha=0.7)
            self.canvas.draw()
            self.status_bar.config(text=f"Referans ayarlandı: ({self.ref_x:.4f}, {self.ref_y:.4f})")
        else:
            messagebox.showinfo("Referans", "Önce grafik üzerinde bir noktaya tıklayın!")
    
    def clear_reference_line(self):
        """Referans çizgilerini temizle"""
        self.ref_x = None
        self.ref_y = None
        self.ref_label.config(text="")
        self.update_plots()
    
    def on_mouse_move(self, event):
        """Fare hareketinde koordinatları göster"""
        if event.inaxes:
            self.coord_label.config(text=f"x: {event.xdata:.4f}, y: {event.ydata:.4f}")
            
            # Referans farkını göster
            if hasattr(self, 'ref_x') and self.ref_x:
                dx = event.xdata - self.ref_x
                dy = event.ydata - self.ref_y
                self.ref_label.config(text=f"Δx: {dx:.4f}, Δy: {dy:.4f}")
    
    def create_file_section(self, parent):
        """Dosya seçimi bölümü"""
        file_frame = ttk.LabelFrame(parent, text="📂 Veri Kaynağı", padding="8")
        file_frame.pack(fill=tk.X, pady=(0, 8))
        
        self.file_path = tk.StringVar()
        ttk.Entry(file_frame, textvariable=self.file_path, width=30).pack(fill=tk.X, pady=2)
        
        btn_frame = ttk.Frame(file_frame)
        btn_frame.pack(fill=tk.X, pady=2)
        
        ttk.Button(btn_frame, text="CSV Aç", command=self.browse_file, width=10).pack(side=tk.LEFT, padx=2)
        ttk.Button(btn_frame, text="PNG Aç", command=self.browse_image, width=10).pack(side=tk.LEFT, padx=2)
        ttk.Button(btn_frame, text="Kalibrasyon", command=self.show_calibration_dialog, width=12).pack(side=tk.LEFT, padx=2)
        
        self.file_status = ttk.Label(file_frame, text="Dosya yok", style='Status.TLabel')
        self.file_status.pack(pady=(5, 0))
    
    def create_channel_section(self, parent, channel, color):
        """Kanal ayarları bölümü"""
        frame = ttk.LabelFrame(parent, text=f"{channel} Ayarları", padding="8")
        frame.pack(fill=tk.X, pady=(0, 8))
        
        color_hex = {'CH1': '#00d4ff', 'CH2': '#ff8844'}[channel]
        
        # Enable checkbox
        var_enabled = tk.BooleanVar(value=(channel == 'CH1'))
        setattr(self, f'{channel.lower()}_enabled', var_enabled)
        ttk.Checkbutton(frame, text=f"{channel} Aktif", variable=var_enabled).pack(anchor='w')
        
        # Tip seçimi
        type_var = tk.StringVar(value="Akim")
        setattr(self, f'{channel.lower()}_type', type_var)
        ttk.Label(frame, text="Tip:").pack(anchor='w')
        ttk.Combobox(frame, textvariable=type_var, values=['Akim', 'Voltaj'], width=12).pack(anchor='w', pady=(0, 5))
        
        # Ratio preset
        preset_var = tk.StringVar(value="5A->0.25V (20 A/V)")
        setattr(self, f'{channel.lower()}_preset', preset_var)
        ttk.Label(frame, text="Ratio Preset:").pack(anchor='w')
        preset_cb = ttk.Combobox(frame, textvariable=preset_var, values=list(RATIO_PRESETS.keys()), width=18)
        preset_cb.pack(anchor='w', pady=(0, 2))
        preset_cb.bind('<<ComboboxSelected>>', lambda e: self.update_ratio(channel))
        
        # Manuel ratio
        ratio_var = tk.StringVar(value="20.0")
        setattr(self, f'{channel.lower()}_ratio', ratio_var)
        ttk.Label(frame, text="Manuel Ratio (A/V):").pack(anchor='w')
        ttk.Entry(frame, textvariable=ratio_var, width=10).pack(anchor='w')
        
        # Filtre ayarları
        filter_frame = ttk.Frame(frame)
        filter_frame.pack(fill=tk.X, pady=(5, 0))
        
        filter_enabled = tk.BooleanVar(value=False)
        setattr(self, f'{channel.lower()}_filter_enabled', filter_enabled)
        ttk.Checkbutton(filter_frame, text="Filtre", variable=filter_enabled).pack(side=tk.LEFT)
        
        filter_type = tk.StringVar(value="savgol")
        setattr(self, f'{channel.lower()}_filter_type', filter_type)
        ttk.Combobox(filter_frame, textvariable=filter_type, values=['lowpass', 'savgol', 'moving_avg'], width=8).pack(side=tk.LEFT, padx=2)
        
        cutoff_var = tk.StringVar(value="2500")
        setattr(self, f'{channel.lower()}_filter_cutoff', cutoff_var)
        ttk.Entry(filter_frame, textvariable=cutoff_var, width=6).pack(side=tk.LEFT, padx=2)
    
    def create_analysis_section(self, parent):
        """Analiz ayarları bölümü"""
        frame = ttk.LabelFrame(parent, text="⚙️ Analiz Ayarları", padding="8")
        frame.pack(fill=tk.X, pady=(0, 8))
        
        # Preset seçimi
        preset_var = tk.StringVar(value="IEC61000-3-2 Class A")
        self.analysis_preset = preset_var
        ttk.Label(frame, text="Analiz Preset:").pack(anchor='w')
        ttk.Combobox(frame, textvariable=preset_var, values=list(ANALYSIS_PRESETS.keys()), width=20).pack(anchor='w')
        
        # Harmonik sayısı
        harm_var = tk.StringVar(value="40")
        self.num_harmonics = harm_var
        ttk.Label(frame, text="Max Harmonik:").pack(anchor='w', pady=(5, 0))
        ttk.Spinbox(frame, from_=10, to=50, textvariable=harm_var, width=10).pack(anchor='w')
        
        # View mode
        ttk.Label(frame, text="Görünüm:").pack(anchor='w', pady=(5, 0))
        self.view_mode = tk.StringVar(value='overlay')
        view_frame = ttk.Frame(frame)
        view_frame.pack(anchor='w')
        ttk.Radiobutton(view_frame, text="Overlay", variable=self.view_mode, value='overlay').pack(side=tk.LEFT)
        ttk.Radiobutton(view_frame, text="Ayrı", variable=self.view_mode, value='separate').pack(side=tk.LEFT)
        ttk.Radiobutton(view_frame, text="Karşılaştır", variable=self.view_mode, value='compare').pack(side=tk.LEFT)

        # CH1-CH2 Fark Grafiği Filtresi
        ttk.Label(frame, text="CH1-CH2 Fark Filtresi:").pack(anchor='w', pady=(8, 0))
        diff_filter_frame = ttk.Frame(frame)
        diff_filter_frame.pack(fill=tk.X)

        self.diff_filter_enabled = tk.BooleanVar(value=False)
        ttk.Checkbutton(diff_filter_frame, text="Filtre", variable=self.diff_filter_enabled).pack(side=tk.LEFT)

        self.diff_filter_type = tk.StringVar(value="savgol")
        ttk.Combobox(diff_filter_frame, textvariable=self.diff_filter_type,
                     values=['lowpass', 'savgol', 'moving_avg'], width=8).pack(side=tk.LEFT, padx=2)

        self.diff_filter_cutoff = tk.StringVar(value="500")
        ttk.Entry(diff_filter_frame, textvariable=self.diff_filter_cutoff, width=6).pack(side=tk.LEFT, padx=2)
    
    def create_action_section(self, parent):
        """İşlem butonları bölümü"""
        frame = ttk.Frame(parent)
        frame.pack(fill=tk.X, pady=(10, 0))
        
        ttk.Button(frame, text="🔬 ANALİZ ET", command=self.run_analysis, width=15).pack(fill=tk.X, pady=2)
        ttk.Button(frame, text="📄 Rapor Kaydet", command=self.save_report, width=15).pack(fill=tk.X, pady=2)
        ttk.Button(frame, text="💾 PNG Kaydet", command=self.save_figure, width=15).pack(fill=tk.X, pady=2)
        ttk.Button(frame, text="📊 CSV Export", command=self.export_csv, width=15).pack(fill=tk.X, pady=2)
    
    def setup_batch_tab(self):
        """Batch işlem sekmesi"""
        control_frame = ttk.Frame(self.batch_tab, style='Card.TFrame', padding="15")
        control_frame.pack(fill=tk.X, pady=10)
        
        ttk.Label(control_frame, text="Toplu Dosya İşleme", style='Title.TLabel').pack(pady=(0, 10))
        
        ttk.Button(control_frame, text="📁 Dosya Ekle", command=self.batch_add_files, width=20).pack(fill=tk.X, pady=5)
        ttk.Button(control_frame, text="▶️ Batch Analiz Başlat", command=self.run_batch_analysis, width=20).pack(fill=tk.X, pady=5)
        ttk.Button(control_frame, text="📄 Tüm Raporu Kaydet", command=self.save_batch_report, width=20).pack(fill=tk.X, pady=5)
        
        # Dosya listesi
        list_frame = ttk.LabelFrame(self.batch_tab, text="Dosya Listesi", padding="10")
        list_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)
        
        self.batch_listbox = tk.Listbox(list_frame, bg='#16213e', fg='white', selectbackground='#00d4ff')
        self.batch_listbox.pack(fill=tk.BOTH, expand=True)
        
        # Progress
        self.batch_progress = ttk.Progressbar(self.batch_tab, mode='determinate')
        self.batch_progress.pack(fill=tk.X, padx=10, pady=10)
        
        self.batch_status = ttk.Label(self.batch_tab, text="Hazır", style='Status.TLabel')
        self.batch_status.pack()
    
    def setup_report_tab(self):
        """Rapor sekmesi"""
        report_frame = ttk.Frame(self.report_tab, style='Dark.TFrame', padding="10")
        report_frame.pack(fill=tk.BOTH, expand=True)
        
        self.report_text = tk.Text(report_frame, font=('Consolas', 10),
                                   bg='#16213e', fg='#e8e8e8', relief='flat')
        self.report_text.pack(fill=tk.BOTH, expand=True)
        
        report_btn_frame = ttk.Frame(report_frame)
        report_btn_frame.pack(fill=tk.X, pady=10)
        
        ttk.Button(report_btn_frame, text="🔄 Yenile", command=self.refresh_report).pack(side=tk.LEFT, padx=5)
        ttk.Button(report_btn_frame, text="📄 Kaydet", command=self.save_report).pack(side=tk.LEFT, padx=5)
        ttk.Button(report_btn_frame, text="📋 Kopyala", command=self.copy_report).pack(side=tk.LEFT, padx=5)
    
    def setup_shortcuts(self):
        """Klavye kısayolları"""
        self.root.bind('<Control-o>', lambda e: self.browse_file())
        self.root.bind('<Control-i>', lambda e: self.browse_image())
        self.root.bind('<F5>', lambda e: self.run_analysis())
        self.root.bind('<F1>', lambda e: self.show_help())
    
    def show_help(self):
        """Yardım dialog"""
        help_text = """
🔬 PROFESYONEL HARMONİK ANALİZÖR - YARDIM

KLAVYE KISAYOLLARI:
  Ctrl+O   - CSV Dosyası Aç
  Ctrl+I   - PNG Görüntü Aç
  F5       - Analiz Et
  F1       - Bu Yardım

KULLANIM:
1. CSV dosyası veya PNG görüntü yükleyin
2. Kanal ayarlarını yapın (ratio, tip)
3. İsteğe göre filtre uygulayın
4. "Analiz Et" butonuna basın
5. Sonuçları grafiklerde ve raporda görün

PNG KALİBRASYON:
  - Grid sınırlarını ayarlayın
  - X ve Y eksen ölçeklerini belirtin
  - Otomatik veya manuel kalibrasyon

EXPORT:
  - PNG: Grafik görüntüsü
  - CSV: Harmonik verileri
  - TXT: Tam rapor
"""
        messagebox.showinfo("Yardım", help_text)
    
    def update_ratio(self, channel):
        """Ratio preset güncelle"""
        preset = getattr(self, f'{channel.lower()}_preset').get()
        value = RATIO_PRESETS.get(preset)
        if value is not None:
            getattr(self, f'{channel.lower()}_ratio').set(str(value))
    
    def browse_file(self):
        """CSV dosyası seç"""
        filepath = filedialog.askopenfilename(
            title="CSV Dosyası Seç",
            filetypes=[("CSV", "*.csv"), ("All", "*.*")],
            initialdir=os.getcwd()
        )
        if filepath:
            self.file_path.set(filepath)
            self.load_file(filepath)
    
    def browse_image(self):
        """PNG görüntü seç"""
        filepath = filedialog.askopenfilename(
            title="PNG Görüntü Seç",
            filetypes=[("PNG", "*.png"), ("JPG", "*.jpg"), ("All", "*.*")],
            initialdir=os.getcwd()
        )
        if filepath:
            self.file_path.set(filepath)
            self.load_image(filepath)
    
    def show_calibration_dialog(self):
        """Kalibrasyon dialog"""
        dialog = tk.Toplevel(self.root)
        dialog.title("PNG Kalibrasyon")
        dialog.geometry("400x500")
        dialog.configure(bg='#1a1a2e')
        
        calib = self.image_extractor.calibration.copy()
        
        # X eksen
        ttk.Label(dialog, text="X Ekseni", style='Title.TLabel').pack(pady=(10, 5))
        ttk.Label(dialog, text="X Başlangıç (piksel):").pack()
        x0_var = tk.StringVar(value=str(calib['x0']))
        ttk.Entry(dialog, textvariable=x0_var).pack()
        ttk.Label(dialog, text="X Bitiş (piksel):").pack()
        x1_var = tk.StringVar(value=str(calib['x1']))
        ttk.Entry(dialog, textvariable=x1_var).pack()
        
        # Y eksen
        ttk.Label(dialog, text="Y Ekseni", style='Title.TLabel').pack(pady=(10, 5))
        ttk.Label(dialog, text="Y Başlangıç (piksel):").pack()
        y0_var = tk.StringVar(value=str(calib['y0']))
        ttk.Entry(dialog, textvariable=y0_var).pack()
        ttk.Label(dialog, text="Y Bitiş (piksel):").pack()
        y1_var = tk.StringVar(value=str(calib['y1']))
        ttk.Entry(dialog, textvariable=y1_var).pack()
        
        # Ölçekler
        ttk.Label(dialog, text="Ölçekler", style='Title.TLabel').pack(pady=(10, 5))
        ttk.Label(dialog, text="Zaman Ölçeği (s):").pack()
        time_var = tk.StringVar(value=str(calib['time_scale']))
        ttk.Entry(dialog, textvariable=time_var).pack()
        
        ttk.Label(dialog, text="Voltaj Ölçeği (V):").pack()
        volt_var = tk.StringVar(value=str(calib['volt_scale']))
        ttk.Entry(dialog, textvariable=volt_var).pack()
        
        def save_calib():
            try:
                self.image_extractor.calibration = {
                    'x0': int(x0_var.get()), 'x1': int(x1_var.get()),
                    'y0': int(y0_var.get()), 'y1': int(y1_var.get()),
                    'time_scale': float(time_var.get()),
                    'volt_scale': float(volt_var.get()),
                    'time_unit': 's', 'volt_unit': 'V'
                }
                dialog.destroy()
                messagebox.showinfo("Başarılı", "Kalibrasyon kaydedildi")
            except ValueError:
                messagebox.showerror("Hata", "Geçersiz değer")
        
        ttk.Button(dialog, text="Kaydet", command=save_calib).pack(pady=20)
    
    def load_file(self, filepath):
        """CSV dosyası yükle"""
        try:
            with open(filepath, 'r') as f:
                # Dosya formatını kontrol et. Dalga formu CSV'leri "Model:" ile başlamaz.
                line1 = f.readline()
                if line1.startswith("Model:"):
                    raise ValueError("Bu bir ayar dosyası gibi görünüyor, dalga formu verisi değil. Lütfen osiloskoptan dalga formunu CSV olarak kaydedin.")
                f.seek(0) # Dosyayı başa sar

                header1_str = f.readline().strip()
                header2_str = f.readline().strip()

                header1 = header1_str.split(',')
                header2 = header2_str.split(',')
            
            has_ch1 = 'CH1' in header1
            has_ch2 = 'CH2' in header1
            print(f"Kanal tespiti: CH1={has_ch1}, CH2={has_ch2}")
            
            try:
                if has_ch1 and has_ch2:
                    print("Çift kanal (CH1 & CH2) modu...")
                    start_time = float(header2[3])
                    increment = float(header2[4])
                    df = pd.read_csv(filepath, skiprows=2, header=None, usecols=[0, 1, 2],
                                    names=['index', 'ch1', 'ch2'])
                    ch1_data = df['ch1'].values
                    ch2_data = df['ch2'].values
                elif has_ch1:
                    print("Tek kanal (CH1) modu...")
                    start_time = float(header2[2])
                    increment = float(header2[3])
                    df = pd.read_csv(filepath, skiprows=2, header=None, usecols=[0, 1],
                                    names=['index', 'ch1'])
                    ch1_data = df['ch1'].values
                    ch2_data = None
                elif has_ch2:
                    print("Tek kanal (CH2) modu...")
                    start_time = float(header2[2])
                    increment = float(header2[3])
                    df = pd.read_csv(filepath, skiprows=2, header=None, usecols=[0, 1],
                                    names=['index', 'ch2'])
                    ch1_data = None
                    ch2_data = df['ch2'].values
                else:
                    raise ValueError("CSV dosyasında CH1 veya CH2 kanalı bulunamadı.")

            except IndexError as ie:
                raise ValueError(f"CSV başlık formatı hatalı (IndexError). Beklenen Rigol dalga formu formatında değil. Header2: '{header2_str}'")
            
            time = start_time + df['index'].values * increment
            sample_rate = 1 / increment
            
            self.data = {
                'time': time,
                'ch1': ch1_data,
                'ch2': ch2_data,
                'dt': increment,
                'sample_rate': sample_rate,
                'has_ch2': has_ch2,
                'filepath': filepath,
                'source': 'csv'
            }
            
            n_points = len(df)
            duration = n_points * increment * 1000
            fname = os.path.basename(filepath)
            
            status = f"✓ {fname} | {n_points:,} nokta | {duration:.1f}ms"
            if has_ch1 and has_ch2: status += " | CH1+CH2"
            elif has_ch1: status += " | CH1"
            elif has_ch2: status += " | CH2"
            self.file_status.config(text=status, foreground="#00ff88")
            
            self.ch1_enabled.set(has_ch1)
            self.ch2_enabled.set(has_ch2)
            
            self.status_bar.config(text=f"Yüklü: {os.path.basename(filepath)}")
            
        except ValueError as e:
            # Kendi oluşturduğumuz veya formatla ilgili ValueError'ları yakala
            self.file_status.config(text=f"Hata: {str(e)}", foreground="#ff4444")
            messagebox.showerror("Dosya Yükleme Hatası", str(e))
        except Exception as e:
            self.file_status.config(text=f"CSV okuma hatası: {str(e)}. Format uyumsuz.", foreground="#ff4444")
            messagebox.showerror("Dosya Yükleme Hatası", f"CSV dosyası okunurken genel bir hata oluştu: {e}\n\nLütfen dalga formu verisi içeren geçerli bir Rigol CSV dosyası seçtiğinizden emin olun.")
    
    def load_image(self, filepath):
        """PNG görüntü yükle ve dalga formu çıkar"""
        result = self.image_extractor.extract_waveform(filepath)
        
        if result is not None:
            # PNG'den çıkarılan veriyi standardize et
            n_samples = len(result['signal'])
            sample_rate = result['sample_rate']
            time = np.linspace(0, n_samples / sample_rate, n_samples)
            
            self.data = {
                'time': time,
                'ch1': result['signal'],
                'ch2': None,
                'dt': 1 / sample_rate,
                'sample_rate': sample_rate,
                'has_ch2': False,
                'filepath': filepath,
                'source': 'png',
                'calibration': result.get('calibration', {})
            }
            
            status = f"✓ PNG: {n_samples:,} nokta, {time[-1]*1000:.1f}ms, {sample_rate/1e3:.1f}kHz"
            self.file_status.config(text=status, foreground="#00d4ff")
            self.ch2_enabled.set(False)
            
            self.status_bar.config(text=f"Yüklü: {os.path.basename(filepath)} (PNG çıkarıldı)")
        else:
            self.file_status.config(text="Görüntü işleme hatası", foreground="#ff4444")
    
    def apply_filter(self, signal, sample_rate, channel):
        """Filtre uygula"""
        if channel == 'CH1':
            if not self.ch1_filter_enabled.get():
                return signal, False, ""
            filter_type = self.ch1_filter_type.get()
            try:
                cutoff = float(self.ch1_filter_cutoff.get())
            except:
                cutoff = 2500
        else:
            if not self.ch2_filter_enabled.get():
                return signal, False, ""
            filter_type = self.ch2_filter_type.get()
            try:
                cutoff = float(self.ch2_filter_cutoff.get())
            except:
                cutoff = 2500
        
        filter_info = f" | {filter_type}"
        
        if filter_type == 'lowpass':
            nyq = sample_rate / 2
            cutoff = min(cutoff, nyq * 0.9)
            b, a = butter(4, cutoff / nyq, btype='low')
            filter_info += f" {cutoff:.0f}Hz"
            return filtfilt(b, a, signal), True, filter_info
        
        elif filter_type == 'savgol':
            window = 51  # Must be an odd number
            filter_info += f" w={window}"
            return savgol_filter(signal, window, 3), True, filter_info
        
        elif filter_type == 'moving_avg':
            window = 51
            kernel = np.ones(window) / window
            filter_info += f" w={window}"
            return np.convolve(signal, kernel, mode='same'), True, filter_info
        
        return signal, False, ""

    def apply_diff_filter(self, signal, sample_rate):
        """CH1-CH2 fark grafiği için filtre uygula"""
        filter_type = self.diff_filter_type.get()
        try:
            cutoff = float(self.diff_filter_cutoff.get())
        except:
            cutoff = 500

        filter_label = f" [{filter_type}"

        if filter_type == 'lowpass':
            nyq = sample_rate / 2
            cutoff = min(cutoff, nyq * 0.9)
            b, a = butter(4, cutoff / nyq, btype='low')
            filter_label += f" {cutoff:.0f}Hz]"
            return filtfilt(b, a, signal), filter_label

        elif filter_type == 'savgol':
            window = int(cutoff) if cutoff > 10 else 51
            if window % 2 == 0:
                window += 1  # Must be odd
            window = min(window, len(signal) - 1)
            if window < 5:
                window = 5
            filter_label += f" w={window}]"
            return savgol_filter(signal, window, 3), filter_label

        elif filter_type == 'moving_avg':
            window = int(cutoff) if cutoff > 1 else 51
            window = min(window, len(signal) - 1)
            kernel = np.ones(window) / window
            filter_label += f" w={window}]"
            return np.convolve(signal, kernel, mode='same'), filter_label

        return signal, ""

    def run_analysis(self):
        """Ana analiz fonksiyonu"""
        if not self.data:
            messagebox.showwarning("Uyarı", "Lütfen veri yükleyin!")
            return
        
        try:
            num_harm = int(self.num_harmonics.get())
        except:
            num_harm = 40
        
        self.results = {}
        time = self.data['time']
        sample_rate = self.data['sample_rate']
        
        # CH1 Analizi
        if self.ch1_enabled.get() and self.data['ch1'] is not None:
            try:
                ratio = float(self.ch1_ratio.get())
            except:
                ratio = 20.0
            
            ch_type = self.ch1_type.get()
            raw_data = self.data['ch1']
            
            if ch_type == 'Akim':
                signal = raw_data * ratio
                unit = 'A'
            else:
                signal = raw_data * 10
                unit = 'V'
            
            signal_filtered, filter_active, filter_info = self.apply_filter(signal, sample_rate, 'CH1')
            
            metrics = self.analyzer.calculate_all_metrics(signal_filtered, sample_rate, num_harmonics=num_harm)
            
            self.results['CH1'] = {
                'channel': 'CH1',
                'type': ch_type,
                'unit': unit,
                'ratio': ratio,
                'time': time[:len(signal_filtered)],
                'signal': signal_filtered,
                'signal_raw': signal,
                'sample_rate': sample_rate,
                'filter_active': filter_active,
                'filter_info': filter_info,
                **metrics
            }
        
        # CH2 Analizi
        if self.ch2_enabled.get() and self.data.get('ch2') is not None:
            try:
                ratio = float(self.ch2_ratio.get())
            except:
                ratio = 20.0
            
            ch_type = self.ch2_type.get()
            raw_data = self.data['ch2']
            
            if ch_type == 'Akim':
                signal = raw_data * ratio
                unit = 'A'
            else:
                signal = raw_data * 10
                unit = 'V'
            
            signal_filtered, filter_active, filter_info = self.apply_filter(signal, sample_rate, 'CH2')
            
            metrics = self.analyzer.calculate_all_metrics(signal_filtered, sample_rate, num_harmonics=num_harm)
            
            self.results['CH2'] = {
                'channel': 'CH2',
                'type': ch_type,
                'unit': unit,
                'ratio': ratio,
                'time': time[:len(signal_filtered)],
                'signal': signal_filtered,
                'signal_raw': signal,
                'sample_rate': sample_rate,
                'filter_active': filter_active,
                'filter_info': filter_info,
                **metrics
            }

        # CH1-CH2 FARK ANALİZİ
        if 'CH1' in self.results and 'CH2' in self.results:
            ch1_res = self.results['CH1']
            ch2_res = self.results['CH2']

            # Fark sinyali oluştur
            min_len = min(len(ch1_res['signal']), len(ch2_res['signal']))
            diff_signal = ch1_res['signal'][:min_len] - ch2_res['signal'][:min_len]
            diff_time = ch1_res['time'][:min_len]

            # Fark sinyaline filtre uygula (opsiyonel)
            filter_info_diff = ''
            if self.diff_filter_enabled.get():
                diff_signal, filter_info_diff = self.apply_diff_filter(diff_signal, sample_rate)

            # Fark sinyalinin tam analizi
            diff_metrics = self.analyzer.calculate_all_metrics(diff_signal, sample_rate, num_harmonics=num_harm)

            # Birim belirleme (her iki kanal aynı türse o birim, değilse genel)
            if ch1_res['type'] == ch2_res['type']:
                diff_unit = ch1_res['unit']
                diff_type = ch1_res['type']
            else:
                diff_unit = 'V/A'
                diff_type = 'Karma'

            self.results['DIFF'] = {
                'channel': 'CH1-CH2',
                'type': diff_type,
                'unit': diff_unit,
                'ratio': 1.0,
                'time': diff_time,
                'signal': diff_signal,
                'signal_raw': diff_signal,
                'sample_rate': sample_rate,
                'filter_active': self.diff_filter_enabled.get(),
                'filter_info': filter_info_diff,
                **diff_metrics
            }

        self.update_plots()
        self.display_results()
        self.refresh_report()
        
        self.status_bar.config(text=f"Analiz tamamlandı | {len(self.results)} kanal | {datetime.now().strftime('%H:%M:%S')}")
    
    def on_canvas_click(self, event):
        """Canvas tıklama - cursor ölçümü"""
        if event.inaxes is None:
            return
        
        if self.active_cursor is None:
            self.cursor1_pos = (event.xdata, event.ydata)
            self.active_cursor = 1
            self.status_bar.config(text=f"Cursor 1: x={event.xdata:.4f}, y={event.ydata:.4f}")
        else:
            self.cursor2_pos = (event.xdata, event.ydata)
            self.active_cursor = None
            
            if self.cursor1_pos and self.cursor2_pos:
                dx = abs(self.cursor2_pos[0] - self.cursor1_pos[0])
                dy = abs(self.cursor2_pos[1] - self.cursor1_pos[1])
                self.status_bar.config(text=f"Δx={dx:.4f}, Δy={dy:.4f}")
            
            self.update_plots()
    
    def update_plots(self):
        """Grafikleri güncelle"""
        if not self.results:
            return
        
        self.fig.clear()
        mode = self.view_mode.get()
        
        if mode == 'separate':
            self.plot_separate()
        elif mode == 'compare':
            self.plot_compare()
        else:
            self.plot_overlay()
        
        self.fig.tight_layout()
        self.canvas.draw()
    
    def plot_overlay(self):
        """Overlay grafik"""
        ax1 = self.fig.add_subplot(3, 2, 1)
        ax2 = self.fig.add_subplot(3, 2, 2)
        ax3 = self.fig.add_subplot(3, 2, 3)
        ax4 = self.fig.add_subplot(3, 2, 4)
        ax5 = self.fig.add_subplot(3, 2, 5)
        ax6 = self.fig.add_subplot(3, 2, 6)

        # DIFF hariç ana kanallar
        main_channels = {ch: res for ch, res in self.results.items() if ch != 'DIFF'}
        colors = {'CH1': '#00d4ff', 'CH2': '#ff8844'}
        width = 0.35
        h_nums = list(range(1, 41))

        # Harmonik bar
        ax1.set_facecolor('#16213e')
        for i, (ch, res) in enumerate(main_channels.items()):
            amps = [h['amplitude'] * 1000 for h in res['harmonics'][:40]]
            offset = -width/2 if i == 0 else width/2
            label = f'{ch} THD={res["thd"]:.1f}%'
            ax1.bar([x + offset for x in h_nums], amps, width, color=colors.get(ch, '#ffffff'), alpha=0.7, label=label)

        any_current = any(res['type'] == 'Akim' for ch, res in main_channels.items())
        if any_current:
            limits = [IEC_CLASS_A_LIMITS.get(h, 0) * 1000 for h in range(1, 41)]
            ax1.step(h_nums[1:], limits[1:], where='mid', color='#ffaa00', linewidth=2, linestyle='--', label='IEC Limit')

        ax1.set_xlabel('Harmonik No', color='white')
        ax1.set_ylabel('Genlik (mA/mV)', color='white')
        ax1.set_title('Harmonik Spektrum', color='white', fontweight='bold')
        ax1.set_xlim(0, 42)
        ax1.legend(loc='upper right', facecolor='#16213e', labelcolor='white', fontsize=8)
        ax1.grid(True, alpha=0.3, axis='y')
        ax1.tick_params(colors='white')

        # Dalga formu 60ms
        ax2.set_facecolor('#16213e')
        for ch, res in main_channels.items():
            samples_60ms = int(0.060 * res['sample_rate'])
            samples_to_show = min(samples_60ms, len(res['signal']))
            start = len(res['signal']) // 2 - samples_to_show // 2
            time_ms = (res['time'][start:start+samples_to_show] - res['time'][start]) * 1000
            ax2.plot(time_ms, res['signal'][start:start+samples_to_show] * 1000,
                    color=colors.get(ch, '#ffffff'), linewidth=0.6, label=f'{ch} RMS={res["rms"]*1000:.1f}m{res["unit"]}')
        ax2.axhline(0, color='gray', linestyle='--', linewidth=0.5)
        ax2.set_xlabel('Zaman (ms)', color='white')
        ax2.set_ylabel('Genlik (mA/mV)', color='white')
        ax2.set_title('Dalga Formu (60ms)', color='white')
        ax2.legend(loc='upper right', facecolor='#16213e', labelcolor='white', fontsize=8)
        ax2.grid(True, alpha=0.3)
        ax2.tick_params(colors='white')
        
        # Limit yüzdeleri
        ax3.set_facecolor('#16213e')
        current_channels = {ch: res for ch, res in main_channels.items() if res['type'] == 'Akim'}
        if current_channels:
            for i, (ch, res) in enumerate(current_channels.items()):
                percents = [h['percent'] for h in res['harmonics'][1:41]]
                offset = -width/2 if i == 0 else width/2
                bar_colors = [colors.get(ch, '#ffffff') if p <= 100 else '#ff4444' for p in percents]
                ax3.bar([x + offset for x in range(2, 41)], percents, width, color=bar_colors, alpha=0.7, label=ch)
            ax3.axhline(100, color='red', linestyle='--', linewidth=2, label='100% Limit')
            ax3.set_xlabel('Harmonik No', color='white')
            ax3.set_ylabel('Limite Göre (%)', color='white')
            ax3.set_title('IEC Limit Karşılaştırma', color='white')
            ax3.set_xlim(1, 41)
            ax3.legend(loc='upper right', facecolor='#16213e', labelcolor='white', fontsize=8)
        ax3.grid(True, alpha=0.3, axis='y')
        ax3.tick_params(colors='white')
        
        # Dalga formu 10ms
        ax4.set_facecolor('#16213e')
        for ch, res in main_channels.items():
            samples_10ms = int(0.010 * res['sample_rate'])
            samples_to_show = min(samples_10ms, len(res['signal']))
            start = len(res['signal']) // 2 - samples_to_show // 2
            time_ms = (res['time'][start:start+samples_to_show] - res['time'][start]) * 1000
            ax4.plot(time_ms, res['signal'][start:start+samples_to_show] * 1000,
                    color=colors.get(ch, '#ffffff'), linewidth=0.8, label=f'{ch} Pk={res["ipk"]*1000:.1f}m{res["unit"]}')
        ax4.axhline(0, color='gray', linestyle='--', linewidth=0.5)
        ax4.set_xlabel('Zaman (ms)', color='white')
        ax4.set_ylabel('Genlik (mA/mV)', color='white')
        ax4.set_title('Dalga Formu (10ms)', color='white')
        ax4.legend(loc='upper right', facecolor='#16213e', labelcolor='white', fontsize=8)
        ax4.grid(True, alpha=0.3)
        ax4.tick_params(colors='white')
        
        # CH1-CH2 Fark Sinyalinin Harmonik Analizi
        ax5.set_facecolor('#16213e')
        if 'DIFF' in self.results:
            diff_res = self.results['DIFF']

            # Harmonikleri çiz
            h_nums = [h['harmonic'] for h in diff_res['harmonics'][:40]]
            h_amps = [h['amplitude'] * 1000 for h in diff_res['harmonics'][:40]]
            ax5.bar(h_nums, h_amps, color='#00ff88', edgecolor='white', linewidth=0.3, alpha=0.8)

            # Başlıkta tüm önemli verileri göster
            title = f"FARK Harmonik | RMS={diff_res['rms']*1000:.2f}m{diff_res['unit']} | THD={diff_res['thd']:.1f}% | CF={diff_res['cf']:.2f}"
            ax5.set_xlabel('Harmonik No', color='white')
            ax5.set_ylabel('Genlik (mA)', color='white')
            ax5.set_title(title, color='#00ff88', fontsize=9)
            ax5.set_xlim(0, 42)
        else:
            ax5.text(0.5, 0.5, 'Fark için 2 kanal gerekli', transform=ax5.transAxes,
                    ha='center', va='center', color='gray', fontsize=12)
        ax5.grid(True, alpha=0.3, axis='y')
        ax5.tick_params(colors='white')

        # Kanal farkı - Dalga Formu (CH1 - CH2)
        ax6.set_facecolor('#16213e')
        if 'DIFF' in self.results:
            diff_res = self.results['DIFF']

            # 60ms göster
            samples_60ms = int(0.060 * diff_res['sample_rate'])
            samples_to_show = min(samples_60ms, len(diff_res['signal']))
            start = len(diff_res['signal']) // 2 - samples_to_show // 2

            time_ms = (diff_res['time'][start:start+samples_to_show] - diff_res['time'][start]) * 1000
            diff_wave = diff_res['signal'][start:start+samples_to_show] * 1000

            filter_label = diff_res.get('filter_info', '')
            ax6.plot(time_ms, diff_wave, color='#00ff88', linewidth=0.6, label=f'CH1-CH2{filter_label}')
            ax6.axhline(0, color='white', linestyle='--', linewidth=0.5)
            ax6.set_xlabel('Zaman (ms)', color='white')
            ax6.set_ylabel('Fark (mA/mV)', color='white')

            # Başlıkta Peak ve f0 göster
            title = f"FARK Dalga | Pk={diff_res['ipk']*1000:.2f}m{diff_res['unit']} | f0={diff_res['fundamental']:.2f}Hz"
            ax6.set_title(title, color='#00ff88', fontsize=9)
            ax6.legend(loc='upper right', facecolor='#16213e', labelcolor='white', fontsize=8)
        else:
            ax6.text(0.5, 0.5, 'Fark için 2 kanal gerekli', transform=ax6.transAxes,
                    ha='center', va='center', color='gray', fontsize=12)
        ax6.grid(True, alpha=0.3)
        ax6.tick_params(colors='white')
    
    def plot_separate(self):
        """Ayrı grafikler"""
        # DIFF hariç kanalları al
        channels = {ch: res for ch, res in self.results.items() if ch != 'DIFF'}
        n = len(channels)
        colors = {'CH1': '#00d4ff', 'CH2': '#ff8844'}

        for i, (ch, res) in enumerate(channels.items()):
            ax1 = self.fig.add_subplot(2, n*2, 1 + i*2)
            ax2 = self.fig.add_subplot(2, n*2, 2 + i*2)
            ax3 = self.fig.add_subplot(2, n*2, 3 + n*2 + i*2)
            ax4 = self.fig.add_subplot(2, n*2, 4 + n*2 + i*2)
            
            color = colors[ch]
            
            # Harmonik bar
            ax1.set_facecolor('#16213e')
            h_nums = range(1, 41)
            amps = [h['amplitude'] * 1000 for h in res['harmonics'][:40]]
            bar_colors = [color if h['status'] != 'FAIL' else '#ff4444' for h in res['harmonics'][:40]]
            ax1.bar(h_nums, amps, color=bar_colors, edgecolor='white', linewidth=0.3)
            ax1.set_xlabel('Harmonik', color='white')
            ax1.set_ylabel('mA/mV', color='white')
            ax1.set_title(f'{ch} THD={res["thd"]:.1f}%', color=color, fontweight='bold')
            ax1.set_xlim(0, 42)
            ax1.grid(True, alpha=0.3)
            ax1.tick_params(colors='white')
            
            # Dalga formu
            ax2.set_facecolor('#16213e')
            samples_per_period = int(res['sample_rate'] / res['fundamental'])
            samples_to_show = min(2 * samples_per_period, len(res['signal']))
            start = len(res['signal']) // 2 - samples_to_show // 2
            time_ms = (res['time'][start:start+samples_to_show] - res['time'][start]) * 1000
            ax2.plot(time_ms, res['signal'][start:start+samples_to_show] * 1000, color=color, linewidth=0.8)
            ax2.axhline(0, color='gray', linestyle='--', linewidth=0.5)
            ax2.set_xlabel('Zaman (ms)', color='white')
            ax2.set_ylabel('mA/mV', color='white')
            ax2.set_title(f'{ch} Dalga Formu', color='white')
            ax2.grid(True, alpha=0.3)
            ax2.tick_params(colors='white')
            
            # FFT spektrum
            ax3.set_facecolor('#16213e')
            n_pts = len(res['signal'])
            yf = np.abs(fft(res['signal']))[:n_pts//2] * 2 / n_pts * 1000
            xf = fftfreq(n_pts, 1/res['sample_rate'])[:n_pts//2]
            mask = xf <= 2500
            ax3.plot(xf[mask], yf[mask], color=color, linewidth=0.5)
            ax3.set_xlabel('Frekans (Hz)', color='white')
            ax3.set_ylabel('mA/mV', color='white')
            ax3.set_title(f'{ch} FFT', color='white')
            ax3.grid(True, alpha=0.3)
            ax3.tick_params(colors='white')
            
            # Limit yüzdesi
            ax4.set_facecolor('#16213e')
            if res['type'] == 'Akim':
                percents = [h['percent'] for h in res['harmonics'][1:41]]
                bar_colors = ['#00ff88' if p <= 100 else '#ff4444' for p in percents]
                ax4.bar(range(2, 41), percents, color=bar_colors, edgecolor='white', linewidth=0.3)
                ax4.axhline(100, color='red', linestyle='--', linewidth=2)
                ax4.set_xlabel('Harmonik', color='white')
                ax4.set_ylabel('Limite (%)', color='white')
                ax4.set_title(f'{ch} IEC %', color='white')
                ax4.set_xlim(1, 41)
            ax4.grid(True, alpha=0.3)
            ax4.tick_params(colors='white')
    
    def plot_compare(self):
        """Karşılaştırmalı grafik"""
        colors = {'CH1': '#00d4ff', 'CH2': '#ff8844', 'DIFF': '#00ff88'}

        # DIFF hariç kanalları al
        channels = {ch: res for ch, res in self.results.items() if ch != 'DIFF'}

        if len(channels) == 1:
            ch = list(channels.keys())[0]
            res = channels[ch]
            ax1 = self.fig.add_subplot(2, 1, 1)
            ax2 = self.fig.add_subplot(2, 1, 2)
            
            h_nums = range(1, 41)
            amps = [h['amplitude'] * 1000 for h in res['harmonics'][:40]]
            bar_colors = [colors[ch] if h['status'] != 'FAIL' else '#ff4444' for h in res['harmonics'][:40]]
            ax1.bar(h_nums, amps, color=bar_colors, edgecolor='white', linewidth=0.3)
            ax1.set_xlabel('Harmonik', color='white')
            ax1.set_ylabel('mA/mV', color='white')
            ax1.set_title(f'{ch} Harmonikler', color='white', fontweight='bold')
            ax1.grid(True, alpha=0.3)
            ax1.tick_params(colors='white')
            
            samples_per_period = int(res['sample_rate'] / res['fundamental'])
            samples_to_show = min(2 * samples_per_period, len(res['signal']))
            start = len(res['signal']) // 2 - samples_to_show // 2
            time_ms = (res['time'][start:start+samples_to_show] - res['time'][start]) * 1000
            ax2.plot(time_ms, res['signal'][start:start+samples_to_show] * 1000, color=colors[ch], linewidth=0.8)
            ax2.set_xlabel('Zaman (ms)', color='white')
            ax2.set_ylabel('mA/mV', color='white')
            ax2.set_title(f'{ch} Dalga Formu', color='white')
            ax2.grid(True, alpha=0.3)
            ax2.tick_params(colors='white')
            return
        
        ax1 = self.fig.add_subplot(2, 2, 1)
        ax2 = self.fig.add_subplot(2, 2, 2)
        ax3 = self.fig.add_subplot(2, 1, 2)
        
        # Ayrı harmonikler
        for ch, res in channels.items():
            ax = ax1 if ch == 'CH1' else ax2
            h_nums = range(1, 41)
            amps = [h['amplitude'] * 1000 for h in res['harmonics'][:40]]
            ax.bar(h_nums, amps, color=colors.get(ch, '#ffffff'), alpha=0.7, label=ch)
            ax.set_xlabel('Harmonik', color='white')
            ax.set_ylabel('mA/mV', color='white')
            ax.set_title(f'{ch} THD={res["thd"]:.1f}%', color=colors.get(ch, '#ffffff'), fontweight='bold')
            ax.legend(loc='upper right', facecolor='#16213e', labelcolor='white')
            ax.grid(True, alpha=0.3)
            ax.tick_params(colors='white')
        
        # Tablo
        ax3.axis('off')
        table_data = []
        headers = ['H#', 'CH1(mA)', 'CH2(mA)', 'Fark', 'Limit', 'CH1%', 'CH2%']
        
        ch1 = self.results.get('CH1')
        ch2 = self.results.get('CH2')
        
        for h_idx in range(20):
            h = h_idx + 1
            row = [str(h)]
            
            if ch1:
                row.append(f'{ch1["harmonics"][h_idx]["amplitude"]*1000:.2f}')
            else:
                row.append('-')
            
            if ch2:
                row.append(f'{ch2["harmonics"][h_idx]["amplitude"]*1000:.2f}')
            else:
                row.append('-')
            
            if ch1 and ch2:
                diff = (ch1['harmonics'][h_idx]['amplitude'] - ch2['harmonics'][h_idx]['amplitude']) * 1000
                row.append(f'{diff:+.2f}')
            else:
                row.append('-')
            
            limit = IEC_CLASS_A_LIMITS.get(h, 0)
            row.append(f'{limit*1000:.1f}' if h > 1 else '-')
            
            if ch1 and h > 1:
                row.append(f'{ch1["harmonics"][h_idx]["percent"]:.1f}%')
            else:
                row.append('-')
            
            if ch2 and h > 1:
                row.append(f'{ch2["harmonics"][h_idx]["percent"]:.1f}%')
            else:
                row.append('-')
            
            table_data.append(row)
        
        table = ax3.table(cellText=table_data, colLabels=headers, loc='center',
                         cellLoc='center', colColours=['#0f3460']*7)
        table.auto_set_font_size(False)
        table.set_fontsize(9)
        table.scale(1.2, 1.5)
        
        for (row, col), cell in table.get_celld().items():
            cell.set_facecolor('#16213e')
            cell.set_text_props(color='white')
            if row == 0:
                cell.set_facecolor('#0f3460')
                cell.set_text_props(fontweight='bold', color='#00d4ff')
    
    def display_results(self):
        """Sonuçları göster"""
        current_channels = [ch for ch, res in self.results.items() if res['type'] == 'Akim' and ch != 'DIFF']

        # DIFF özet bilgisi
        diff_info = ""
        if 'DIFF' in self.results:
            diff = self.results['DIFF']
            diff_info = f" | FARK: RMS={diff['rms']*1000:.2f}m{diff['unit']} THD={diff['thd']:.1f}%"

        if current_channels:
            all_passed = all(self.results[ch]['passed'] for ch in current_channels)

            if all_passed:
                self.result_banner.config(text=f"✓ IEC 61000-3-2: PASSED{diff_info}",
                                          foreground='#00ff88', font=('Consolas', 14, 'bold'))
            else:
                failed_channels = [ch for ch in current_channels if not self.results[ch]['passed']]
                self.result_banner.config(text=f"✗ IEC: {' '.join(failed_channels)} FAILED{diff_info}",
                                          foreground='#ff4444', font=('Consolas', 14, 'bold'))
        else:
            self.result_banner.config(text=f"Analiz tamamlandı{diff_info}",
                                      foreground='#00d4ff', font=('Consolas', 12))
    
    def generate_report(self):
        """Rapor oluştur"""
        report = f"""
================================================================================
          PROFESYONEL HARMONİK ANALİZ RAPORU
          IEC 61000-3-2 Uyumluluk Raporu
================================================================================
Tarih: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}
Kaynak: {self.data.get('filepath', 'N/A')}
Analiz: IEC 61000-3-2 Class A Limitleri

"""
        
        for ch, res in self.results.items():
            status = "PASSED" if res['passed'] else "FAILED"
            
            report += f"""
--------------------------------------------------------------------------------
{ch} KANAL ANALİZ SONUÇLARI
--------------------------------------------------------------------------------
Tip: {res['type']}
Dönüşüm Oranı: {res['ratio']} A/V
Temel Frekans: {res['fundamental']:.2f} Hz
RMS Değeri: {res['rms']*1000:.2f} m{res['unit']}
Peak Değeri: {res['ipk']*1000:.2f} m{res['unit']}
Crest Factor: {res['cf']:.3f}
THD: {res['thd']:.2f} %
TDD: {res['tdd']:.2f} %
Power Factor: {res['pf']:.4f}
IEC 61000-3-2: {status}
"""
            
            if res['failed']:
                failed_list = ', '.join([f'H{h["harmonic"]}(%{h["percent"]:.1f})' for h in res['failed']])
                report += f"Limit Aşan Harmonikler: {failed_list}\n"
            
            report += """
Sıra   Frekans(Hz)   Genlik(mA)   Limit(mA)    %%      Faz(°)   Durum
"""
            report += "-" * 75 + "\n"
            
            for h in res['harmonics']:
                if h['harmonic'] == 1:
                    report += f"  {h['harmonic']:2d}   {h['frequency']:8.1f}   {h['amplitude']*1000:10.2f}   "
                    report += f"   ---       ---     {h['phase']:6.1f}    FUND\n"
                else:
                    status_icon = "✓" if h['status'] == 'PASS' else "✗"
                    report += f"  {h['harmonic']:2d}   {h['frequency']:8.1f}   {h['amplitude']*1000:10.2f}   "
                    report += f"{h['limit']*1000:8.2f}   {h['percent']:6.1f}%  {h['phase']:6.1f}   {status_icon}\n"
        
        report += """
================================================================================
                          SONUÇLAR
================================================================================
"""
        
        # Özet
        for ch, res in self.results.items():
            status = "PASS" if res['passed'] else "FAIL"
            report += f"{ch}: THD={res['thd']:.2f}%, TDD={res['tdd']:.2f}%, PF={res['pf']:.4f}, IEC={status}\n"
        
        return report
    
    def refresh_report(self):
        """Rapor sekmesini güncelle"""
        if not self.results:
            self.report_text.delete(1.0, tk.END)
            self.report_text.insert(tk.END, "Henüz analiz yapılmadı.\nVeri yükleyip analiz edin.")
            return
        
        report = self.generate_report()
        self.report_text.delete(1.0, tk.END)
        self.report_text.insert(tk.END, report)
    
    def copy_report(self):
        """Raporu panoya kopyala"""
        report = self.generate_report()
        self.root.clipboard_clear()
        self.root.clipboard_append(report)
        messagebox.showinfo("Kopyalandı", "Rapor panoya kopyalandı")
    
    def save_report(self):
        """Rapor kaydet"""
        if not self.results:
            messagebox.showwarning("Uyarı", "Önce analiz yapın!")
            return
        
        filepath = filedialog.asksaveasfilename(
            defaultextension=".txt",
            filetypes=[("Text", "*.txt")],
            initialfile=f"Harmonic_Report_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
        )
        
        if filepath:
            report = self.generate_report()
            with open(filepath, 'w', encoding='utf-8') as f:
                f.write(report)
            messagebox.showinfo("Başarılı", f"Rapor kaydedildi:\n{filepath}")
    
    def save_figure(self):
        """Grafik kaydet"""
        if not self.results:
            messagebox.showwarning("Uyarı", "Önce analiz yapın!")
            return
        
        filepath = filedialog.asksaveasfilename(
            defaultextension=".png",
            filetypes=[("PNG", "*.png"), ("PDF", "*.pdf")],
            initialfile="harmonic_analysis"
        )
        
        if filepath:
            self.fig.savefig(filepath, dpi=150, bbox_inches='tight', facecolor='#1a1a2e')
            messagebox.showinfo("Başarılı", f"Grafik kaydedildi:\n{filepath}")
    
    def export_csv(self):
        """CSV export"""
        if not self.results:
            messagebox.showwarning("Uyarı", "Önce analiz yapın!")
            return
        
        filepath = filedialog.asksaveasfilename(
            defaultextension=".csv",
            filetypes=[("CSV", "*.csv")],
            initialfile=f"harmonic_data_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
        )
        
        if filepath:
            data = []
            for ch, res in self.results.items():
                for h in res['harmonics']:
                    data.append({
                        'Kanal': ch,
                        'Harmonik': h['harmonic'],
                        'Frekans(Hz)': h['frequency'],
                        'Genlik(mA)': h['amplitude'] * 1000,
                        'Limit(mA)': h['limit'] * 1000,
                        'Limit%': h['percent'],
                        'Faz(°)': h['phase'],
                        'Durum': h['status']
                    })
            
            df = pd.DataFrame(data)
            df.to_csv(filepath, index=False, encoding='utf-8')
            messagebox.showinfo("Başarılı", f"Veriler kaydedildi:\n{filepath}")
    
    # ===================== BATCH PROCESSING =====================
    
    def batch_add_files(self):
        """Batch işlem için dosya ekle"""
        filepaths = filedialog.askopenfilenames(
            title="CSV Dosyaları Seç",
            filetypes=[("CSV", "*.csv"), ("All", "*.*")],
            initialdir=os.getcwd()
        )
        
        for fp in filepaths:
            if fp not in self.batch_files:
                self.batch_files.append(fp)
                self.batch_listbox.insert(tk.END, os.path.basename(fp))
    
    def run_batch_analysis(self):
        """Batch analiz çalıştır"""
        if not self.batch_files:
            messagebox.showwarning("Uyarı", "Önce dosya ekleyin!")
            return
        
        self.batch_index = 0
        self.batch_results = []
        self.batch_progress['maximum'] = len(self.batch_files)
        
        def process_next():
            if self.batch_index < len(self.batch_files):
                fp = self.batch_files[self.batch_index]
                self.batch_status.config(text=f"İşleniyor: {os.path.basename(fp)} ({self.batch_index+1}/{len(self.batch_files)})")
                self.root.update()
                
                # Analiz yap
                self.file_path.set(fp)
                self.load_file(fp)
                self.run_analysis()
                
                # Sonuç kaydet
                self.batch_results.append({
                    'file': fp,
                    'results': self.results.copy()
                })
                
                self.batch_index += 1
                self.batch_progress['value'] = self.batch_index
                
                # Sonraki dosya
                self.root.after(100, process_next)
            else:
                self.batch_status.config(text=f"Tamamlandı! {len(self.batch_files)} dosya işlendi.")
                messagebox.showinfo("Tamamlandı", f"Batch işlem tamamlandı.\n{len(self.batch_files)} dosya işlendi.")
        
        process_next()
    
    def save_batch_report(self):
        """Tüm batch sonuçlarını kaydet"""
        if not hasattr(self, 'batch_results') or not self.batch_results:
            messagebox.showwarning("Uyarı", "Önce batch analiz yapın!")
            return
        
        filepath = filedialog.asksaveasfilename(
            defaultextension=".txt",
            filetypes=[("Text", "*.txt")],
            initialfile=f"Batch_Report_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
        )
        
        if filepath:
            report = "=" * 80 + "\n"
            report += "          TOPLU HARMONİK ANALİZ RAPORU\n"
            report += "          Batch Processing Results\n"
            report += "=" * 80 + "\n\n"
            
            for i, br in enumerate(self.batch_results):
                report += f"\n{'='*60}\n"
                report += f"Dosya {i+1}: {br['file']}\n"
                report += f"{'='*60}\n\n"
                
                for ch, res in br['results'].items():
                    status = "PASS" if res['passed'] else "FAIL"
                    report += f"{ch}: THD={res['thd']:.2f}%, TDD={res['tdd']:.2f}%, PF={res['pf']:.4f}, IEC={status}\n"
                    
                    if res['failed']:
                        report += f"  Limit Aşan: {', '.join([f'H{h['harmonic']}' for h in res['failed']])}\n"
            
            with open(filepath, 'w', encoding='utf-8') as f:
                f.write(report)
            
            messagebox.showinfo("Başarılı", f"Batch rapor kaydedildi:\n{filepath}")


def main():
    root = tk.Tk()
    app = DualCurrentAnalyzer(root)
    root.mainloop()


if __name__ == "__main__":
    main()

    root = tk.Tk()
    app = DualCurrentAnalyzer(root)
    root.mainloop()

if __name__ == "__main__":
    main()
