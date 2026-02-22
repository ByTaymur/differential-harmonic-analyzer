# ⚡ İzole Harmonik Analiz Laboratuvarı: Diferansiyel Akım Yöntemi (KCL) ve Kendin Yap (DIY) LC Filtre ile Şebeke Bağımsız THD Ölçümü

![License](https://img.shields.io/badge/license-MIT-blue.svg)
![Python](https://img.shields.io/badge/python-3.8%2B-blue)
![Hardware](https://img.shields.io/badge/hardware-DIY_LC_Filter-orange)

[cite_start]Bu proje, şebeke gerilimindeki (220V AC) mevcut kirlilikten ve arka plan gürültüsünden etkilenmeden, Non-Linear (Doğrusal Olmayan) yüklerin test edilen cihaz (DUT) şebekeye bastığı saf harmonik emisyonlarını ölçmek amacıyla geliştirilmiş donanım ve yazılım mimarisini içermektedir[cite: 7, 8].

Profesyonel laboratuvarlardaki on binlerce dolarlık "AC Grid Simulator" (Şebeke Simülatörü) cihazlarına alternatif olarak geliştirilen bu sistem; temel fizik yasalarını (Empedans ve Kirchhoff Akım Yasası), yaratıcı donanım hack'lerini (damacana ile sıvı soğutma) ve dijital sinyal işlemeyi (DSP) bir araya getiren tam teşekküllü bir Ar-Ge çalışmasıdır.

---

## 📋 İçindekiler
1. [Projenin Amacı ve Karşılaşılan Temel Problem](#1-projenin-amacı-ve-karşılaşılan-temel-problem)
2. [Teorik Altyapı ve Matematiksel Modeller](#2-teorik-altyapı-ve-matematiksel-modeller)
3. [Empedans Duvarı ve Diferansiyel Akım Yöntemi (KCL)](#3-empedans-duvarı-ve-diferansiyel-akım-yöntemi-kcl)
4. [Mühendislik Çözümleri ve DIY Donanım Hack'leri](#4-mühendislik-çözümleri-ve-diy-donanım-hackleri)
5. [Yazılım Mimarisi (DSP) ve Standartlar](#5-yazılım-mimarisi-dsp-ve-standartlar)
6. [Kurulum ve Kullanım](#6-kurulum-ve-kullanım)

---

## 1. Projenin Amacı ve Karşılaşılan Temel Problem

[cite_start]**Problem (Sebep):** Şebekeden (Grid) hiçbir akım çekilmese dahi, hatta bulunan diğer güç elektroniği cihazları nedeniyle şebeke gerilimi saf bir sinüs dalgası olmaktan uzaktır[cite: 2, 7]. [cite_start]Bu durum, herhangi bir yük olmadan şebekeden alınan örnekte (Şekil 1.1) net bir şekilde görülmektedir[cite: 2].



[cite_start]**Sonuç:** Bir cihazın IEC 61000-3-2 standartlarına uygunluğunu test etmek istediğimizde, şebekenin kendi kirliliği osiloskop ölçümlerini manipüle ederek yanıltıcı sonuçlar verir[cite: 7, 35].

[cite_start]**Çözüm İhtiyacı:** Cihazın şebekeye ne kadar harmonik bastığını kesin olarak ölçebilmek için, öncelikle cihaza **temiz (izole) bir referans gerilimi** sağlanması gerekmektedir[cite: 8].

---

## 2. Teorik Altyapı ve Matematiksel Modeller

### 2.1. LC Alçak Geçiren Filtre (Low Pass Filter) Tasarımı
[cite_start]Şebekedeki yüksek frekanslı gürültüleri engellemek amacıyla bir LC Alçak Geçiren Filtre tasarlanmıştır[cite: 10]. 



[Image of LC low pass filter circuit diagram]


Tasarım parametreleri ve kullanılan malzemeler:
* [cite_start]**Endüktans ($L$):** $45\text{ mH}$ ($45 \times 10^{-3}\text{ H}$) [cite: 5, 12]
* [cite_start]**Kapasitans ($C$):** $44\text{ \mu F}$ (2 adet $22\text{ \mu F}$ kapasitörün paralel bağlanmasıyla elde edilmiştir) [cite: 5, 12, 17]

**1. Kesim Frekansı ($f_c$) Hesabı:**
Filtrenin hangi frekanstan sonrasını engellemeye başlayacağını belirleyen formül:
$$f_c = \frac{1}{2\pi\sqrt{LC}}$$
[cite_start]$$f_c = \frac{1}{2 \cdot 3.1415 \cdot \sqrt{0.045 \cdot 0.000044}} \approx \mathbf{113.1\text{ Hz}}$$ [cite: 12]

[cite_start]*Açıklama:* Bu filtre, 113 Hz üzerindeki gürültüleri sönümleyerek tıpkı bir subwoofer (bas) filtresi mantığıyla çalışır[cite: 12, 13]. [cite_start]Şebekenin $50\text{ Hz}$ temel frekansını geçirirken, harmonikleri bloke ederek şebeke kirliliğinden arındırılmış bir ortam oluşturur[cite: 13].



**2. Karakteristik Empedans ($Z_0$) Hesabı:**
[cite_start]Filtrenin ideal çalışması için devrenin karakteristik empedansı da hesaplanmıştır[cite: 12]:
[cite_start]$$Z_0 = \sqrt{\frac{L}{C}} = \sqrt{\frac{0.045}{0.000044}} \approx \mathbf{32\text{ Ohm}}$$ [cite: 12]

---

## 3. Empedans Duvarı ve Diferansiyel Akım Yöntemi (KCL)

[cite_start]Filtreleme şebeke kirliliğini engellemede başarılı olmuştur; ancak sisteme seri giren yüksek endüktans ($45\text{ mH}$), cihazın ürettiği harmoniklerin şebekeye doğru akışını engelleyen bir reaktans bariyeri oluşturmuştur[cite: 4]. [cite_start]Bu durum, düzgün bir harmonik analizi yapmayı zorlaştırmaktadır[cite: 4].

### 3.1. Frekansa Bağlı Reaktans Analizi ($X_L$ ve $X_C$)
Endüktif ve Kapasitif Reaktans formülleri:
$$X_L = 2 \cdot \pi \cdot f \cdot L$$
$$X_C = \frac{1}{2 \cdot \pi \cdot f \cdot C}$$

1. **$50\text{ Hz}$ Temel Frekans için:**
   * $X_L = 2 \cdot 3.14 \cdot 50 \cdot 0.045 = \mathbf{14.1\text{ }\Omega}$
   * $X_C = \frac{1}{2 \cdot 3.14 \cdot 50 \cdot 0.000044} = \mathbf{72.3\text{ }\Omega}$
   * *Durum:* Şebekeden gelen $50\text{ Hz}$ enerjinin empedansı bobin üzerinde düşüktür, cihaz rahatça beslenir.

2. **$250\text{ Hz}$ (5. Harmonik) için (Cihazın ürettiği gürültü):**
   * $X_L = 2 \cdot 3.14 \cdot 250 \cdot 0.045 = \mathbf{70.7\text{ }\Omega}$
   * $X_C = \frac{1}{2 \cdot 3.14 \cdot 250 \cdot 0.000044} = \mathbf{14.4\text{ }\Omega}$
   * *Durum:* Cihazın ürettiği $250\text{ Hz}$'lik akım şebekeye geri dönmek istediğinde $70.7\text{ }\Omega$ gibi yüksek bir duvarla karşılaşır. Akım en düşük dirençli yolu seçeceği için zorunlu olarak empedansı $14.4\text{ }\Omega$ olan kapasitör hattına yönelir.

### 3.2. Çözüm: Diferansiyel Akım ve Harmonik Döngüsü (KCL)
[cite_start]Harmonik akımların empedans sebebiyle şebekeye gidemeyip kapasitöre saptığı öngörülerek "Fark Akımı" yöntemi uygulanmıştır[cite: 4, 21].



Yukarıdaki devrede görüldüğü üzere, sistemde iki adet akım ölçüm noktası kurgulanmıştır:
* **PROB 2 (A2):** Yük akımını ($I_{Load}$) ölçer.
* **PROB 3 (A3):** Filtre/Kapasitör akımını ($I_{Filtre}$) ölçer.
* **Harmonik Döngüsü:** Yükün ürettiği yüksek frekanslı harmonikler, $I_{Harmonik}$ (Döngü) yolunu izleyerek şebekeye ($L$ bobinine) gitmek yerine kapasitör ($C$) üzerinden devresini tamamlar.

Düğüm (Node) noktasındaki formül:
$$I_{Giris} = I_{Kapasitor} + I_{Cihaz}$$

[cite_start]Cihazın harmonik imzasını bulmak için kullanılan diferansiyel denklem[cite: 22]:
$$I_{Cihaz} = I_{Giris} - I_{Kapasitor}$$

---

## 4. Mühendislik Çözümleri ve DIY Donanım Hack'leri

[cite_start]Proje, düşük maliyetli ve erişilebilir malzemelerin mühendislik pratikleriyle birleştirilmesiyle kurulmuştur[cite: 15, 26]. Laboratuvar ortamı olmadan ev/atölye şartlarında geliştirilen çözümler şunlardır:

### 4.1. Reaktör Olarak Standart Rulo Kablo Kullanımı
* **Tasarım:** Endüstriyel bir reaktör satın almak yerine, piyasada kolayca bulunabilen **$300\text{ metre}$ uzunluğunda, $0.75\text{ mm}$ kesitli bakır kablo** kullanılmıştır[cite: 16]. 
* **Avantaj:** Kablo, makarasından sağılmadan kendi sarmal yapısıyla devrede bırakılarak devasa bir hava nüveli bobin elde edilmiş ve hedeflenen $45\text{ mH}$ değerine bu sayede ulaşılmıştır.

### 4.2. "Damacana" ile Pasif Sıvı Soğutma (Thermal Hack)
* **Problem:** Kullanılan bobinin (300m kablo) üzerinden akım geçtiğinde ciddi bir ısınma problemi ortaya çıkar[cite: 18].
* [cite_start]**Çözüm:** Bobinin ısınmasını engellemek için olağanüstü bir pasif soğutma yöntemi geliştirilmiştir: Bobin, su dolu kesilmiş bir damacana içerisine yerleştirilmiştir[cite: 18]. 



### 4.3. Akım-Gerilim Dönüşümü (CT ve Yük Direnci)
* [cite_start]**Ölçüm Donanımı:** Sisteme 2 adet $5\text{A}$ akım trafosu (CT) entegre edilmiştir[cite: 19]. 
* **Çalışma Mantığı:** Osiloskoplar doğrudan akım okuyamadığı için, akım trafolarının çıkışına yük direnci (Burden Resistor) bağlanmıştır. Akım probu dönüşüm oranı yazılımda **$5\text{A} \rightarrow 0.25\text{V}$ ($20\text{ A/V}$)** olarak tanımlanmıştır[cite: 39].

---

## 5. Yazılım Mimarisi (DSP) ve Standartlar

Osiloskop üzerinden gerilim dalga formu olarak `.csv` formatında kaydedilen diferansiyel sensör verileri, Python dilinde geliştirilen **"Profesyonel Harmonik Analizör - Dual Channel Analyzer"** arayüzü ile işlenmektedir[cite: 23, 28].



**Yazılımın Dijital Sinyal İşleme (DSP) Özellikleri:**
* [cite_start]**Dual Channel Okuma:** Tek `.csv` dosyasından CH1 ve CH2 verilerinin senkronize olarak alınması[cite: 30].
* [cite_start]**DC Offset ve Filtreleme:** Sinyalden DC bileşenin çıkarılması (`signal - np.mean(signal)`) ve Savitzky-Golay / Lowpass gibi dijital filtreleme opsiyonları[cite: 176, 915].
* **Matematiksel Ayrıştırma:** Her iki kanalın (Giriş ve Kapasitör) genlikleri dönüştürüldükten sonra $I_{CH1} - I_{CH2}$ fark sinyalinin yazılımsal olarak hesaplanması[cite: 23, 1031].
* [cite_start]**FFT (Hızlı Fourier Dönüşümü):** Zaman domenindeki sinyalin `scipy.fft` kütüphanesi ile frekans domenine aktarılarak ($45-65\text{ Hz}$ arası temel frekansın bulunması dahil) $40.$ harmoniğe kadar olan spektrumun çıkarılması[cite: 48, 174, 213, 231].
* [cite_start]**Güç Kalitesi Metrikleri:** Sinyal üzerinden bağımsız olarak THD (Total Harmonic Distortion), TDD (Total Demand Distortion), RMS, Crest Factor (ipk/rms) ve Power Factor hesaplamaları yapılmaktadır[cite: 183-194].
* **Otomatik IEC Uyumluluk Testi:** Yazılım içerisinde **IEC 61000-3-2 Class A** standart limitleri tanımlıdır[cite: 56, 57-65]. Analiz edilen cihazın her bir harmoniği limitlerle karşılaştırılarak otomatik **PASS / FAIL** raporu oluşturulur[cite: 35, 252].

---

## 6. Kurulum ve Kullanım

### Gereksinimler
Bu projeyi çalıştırmak için aşağıdaki kütüphanelere sahip Python 3.8 veya daha üstü bir sürüm gereklidir:
* [cite_start]`pandas`, `matplotlib`, `numpy`, `scipy`, `tkinter` (GUI için), ve PNG dalga formu çıkarma özellikleri için `Pillow`, `opencv-python`[cite: 41-50, 112-113].

```bash
pip install pandas numpy scipy matplotlib opencv-python pillow
