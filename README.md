# ⚡ İzole Harmonik Analiz Laboratuvarı: Diferansiyel Akım Yöntemi (KCL) ve Kendin Yap (DIY) LC Filtre ile Şebeke Bağımsız THD Ölçümü

![License](https://img.shields.io/badge/license-MIT-blue.svg)
![Python](https://img.shields.io/badge/python-3.8%2B-blue)
![Hardware](https://img.shields.io/badge/hardware-DIY_LC_Filter-orange)

Bu proje, şebeke gerilimindeki (220V AC) mevcut kirlilikten ve arka plan gürültüsünden etkilenmeden, Non-Linear (Doğrusal Olmayan) yüklerin (test edilen cihaz - DUT) şebekeye bastığı saf harmonik emisyonlarını ölçmek amacıyla geliştirilmiş donanım ve yazılım mimarisini içermektedir.

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

**Problem (Sebep):** Şebekeden (Grid) hiçbir akım çekilmese dahi, hatta bulunan diğer güç elektroniği cihazları nedeniyle şebeke gerilimi saf bir sinüs dalgası olmaktan uzaktır. Bu durum, herhangi bir yük olmadan şebekeden alınan örnekte net bir şekilde görülmektedir.

**Sonuç:** Bir cihazın IEC 61000-3-2 standartlarına uygunluğunu test etmek istediğimizde, şebekenin kendi kirliliği osiloskop ölçümlerini manipüle ederek yanıltıcı sonuçlar verir.

**Çözüm İhtiyacı:** Cihazın şebekeye ne kadar harmonik bastığını kesin olarak ölçebilmek için, öncelikle cihaza **temiz (izole) bir referans gerilimi** sağlanması gerekmektedir.

---

## 2. Teorik Altyapı ve Matematiksel Modeller

### 2.1. LC Alçak Geçiren Filtre (Low Pass Filter) Tasarımı
Şebekedeki yüksek frekanslı gürültüleri engellemek amacıyla bir LC Alçak Geçiren Filtre tasarlanmıştır. 

Tasarım parametreleri ve kullanılan malzemeler:
* **Endüktans (L):** **45 mH** (0.045 H)
* **Kapasitans (C):** **44 µF** (2 adet 22 µF kapasitörün paralel bağlanmasıyla elde edilmiştir)

**1. Kesim Frekansı Hesabı:**
Filtrenin hangi frekanstan sonrasını engellemeye başlayacağını belirleyen formül:
$$f_c = \frac{1}{2\pi\sqrt{LC}}$$

Değerler yerine konulduğunda:
$$f_c = \frac{1}{2 \cdot 3.1415 \cdot \sqrt{0.045 \cdot 0.000044}}$$
**Sonuç:** Yaklaşık **113.1 Hz**

*Açıklama:* Bu filtre, 113 Hz üzerindeki gürültüleri sönümleyerek tıpkı bir subwoofer (bas) filtresi mantığıyla çalışır. Şebekenin 50 Hz temel frekansını geçirirken, harmonikleri bloke ederek şebeke kirliliğinden arındırılmış bir ortam oluşturur.

**2. Karakteristik Empedans Hesabı:**
Filtrenin ideal çalışması için devrenin karakteristik empedansı da hesaplanmıştır:
$$Z_0 = \sqrt{\frac{L}{C}} = \sqrt{\frac{0.045}{0.000044}}$$
**Sonuç:** Yaklaşık **32 Ω**

---

## 3. Empedans Duvarı ve Diferansiyel Akım Yöntemi (KCL)

Filtreleme şebeke kirliliğini engellemede başarılı olmuştur; ancak sisteme seri giren yüksek endüktans (45 mH), cihazın ürettiği harmoniklerin şebekeye doğru akışını engelleyen bir reaktans bariyeri oluşturmuştur. Bu durum, düzgün bir harmonik analizi yapmayı zorlaştırmaktadır.

### 3.1. Frekansa Bağlı Reaktans Analizi
Endüktif ($X_L$) ve Kapasitif ($X_C$) Reaktans formülleri:
$$X_L = 2 \cdot \pi \cdot f \cdot L$$
$$X_C = \frac{1}{2 \cdot \pi \cdot f \cdot C}$$

**1. 50 Hz Temel Frekans için:**
* $X_L$ = 2 · 3.14 · 50 · 0.045 = **14.1 Ω**
* $X_C$ = 1 / (2 · 3.14 · 50 · 0.000044) = **72.3 Ω**
* *Durum:* Şebekeden gelen 50 Hz enerjinin empedansı bobin üzerinde düşüktür, cihaz rahatça beslenir.

**2. 250 Hz (5. Harmonik) için (Cihazın ürettiği gürültü):**
* $X_L$ = 2 · 3.14 · 250 · 0.045 = **70.7 Ω**
* $X_C$ = 1 / (2 · 3.14 · 250 · 0.000044) = **14.4 Ω**
* *Durum:* Cihazın ürettiği 250 Hz'lik akım şebekeye geri dönmek istediğinde **70.7 Ω** gibi yüksek bir duvarla karşılaşır. Akım en düşük dirençli yolu seçeceği için zorunlu olarak empedansı **14.4 Ω** olan kapasitör hattına yönelir.

### 3.2. Çözüm: Diferansiyel Akım ve Harmonik Döngüsü (KCL)
Harmonik akımların empedans sebebiyle şebekeye gidemeyip kapasitöre saptığı öngörülerek "Fark Akımı" yöntemi uygulanmıştır.

Düğüm (Node) noktasındaki formül:
$$I_{Giris} = I_{Kapasitor} + I_{Cihaz}$$

Cihazın harmonik imzasını bulmak için kullanılan diferansiyel denklem:
**Cihaz Harmonikleri = Giriş Akımı - Kapasitör Akımı**

---

## 4. Mühendislik Çözümleri ve DIY Donanım Hack'leri

Proje, düşük maliyetli ve erişilebilir malzemelerin mühendislik pratikleriyle birleştirilmesiyle kurulmuştur. Laboratuvar ortamı olmadan ev/atölye şartlarında geliştirilen çözümler şunlardır:

### 4.1. Reaktör Olarak Standart Rulo Kablo Kullanımı
* **Tasarım:** Endüstriyel bir reaktör satın almak yerine, piyasada kolayca bulunabilen **300 metre uzunluğunda, 0.75 mm kesitli bakır kablo** kullanılmıştır. 
* **Avantaj:** Kablo, makarasından sağılmadan kendi sarmal yapısıyla devrede bırakılarak devasa bir hava nüveli bobin elde edilmiş ve hedeflenen **45 mH** değerine bu sayede ulaşılmıştır.

### 4.2. "Damacana" ile Pasif Sıvı Soğutma (Thermal Hack)
* **Problem:** Kullanılan bobinin (300m kablo) üzerinden akım geçtiğinde ciddi bir ısınma problemi ortaya çıkar.
* **Çözüm:** Bobinin ısınmasını engellemek için olağanüstü bir pasif soğutma yöntemi geliştirilmiştir: Bobin (orijinal yalıtkan poşeti içindeyken), su dolu kesilmiş bir damacana içerisine yerleştirilmiştir. Bu sayede suyun termal kapasitesinden faydalanılarak soğutma sağlanmıştır.

### 4.3. Akım-Gerilim Dönüşümü (CT ve Yük Direnci)
* **Ölçüm Donanımı:** Sisteme 2 adet **5A akım trafosu (CT)** entegre edilmiştir. 
* **Çalışma Mantığı:** Osiloskoplar doğrudan akım okuyamadığı için, akım trafolarının çıkışına **100 Ω (100R)** yük direnci (Burden Resistor) bağlanmıştır. Akım probu dönüşüm oranı yazılımda **5A -> 0.25V (20 A/V)** olarak tanımlanmıştır.

---

## 5. Yazılım Mimarisi (DSP) ve Standartlar

Osiloskop üzerinden gerilim dalga formu olarak `.csv` formatında kaydedilen diferansiyel sensör verileri, Python dilinde geliştirilen **"Profesyonel Harmonik Analizör - Dual Channel Analyzer"** arayüzü ile işlenmektedir.

**Yazılımın Dijital Sinyal İşleme (DSP) Özellikleri:**
* **Dual Channel Okuma:** Tek `.csv` dosyasından CH1 ve CH2 verilerinin senkronize olarak alınması.
* **DC Offset ve Filtreleme:** Sinyalden DC bileşenin çıkarılması ve Savitzky-Golay / Lowpass gibi dijital filtreleme opsiyonları.
* **Matematiksel Ayrıştırma:** Her iki kanalın (Giriş ve Kapasitör) genlikleri dönüştürüldükten sonra "CH1 - CH2" fark sinyalinin yazılımsal olarak hesaplanması.
* **FFT (Hızlı Fourier Dönüşümü):** Zaman domenindeki sinyalin `scipy.fft` kütüphanesi ile frekans domenine aktarılarak 40. harmoniğe kadar olan spektrumun çıkarılması.
* **Güç Kalitesi Metrikleri:** Sinyal üzerinden bağımsız olarak THD (Total Harmonic Distortion), TDD (Total Demand Distortion), RMS, Crest Factor (ipk/rms) ve Power Factor hesaplamaları yapılmaktadır.
* **Otomatik IEC Uyumluluk Testi:** Yazılım içerisinde **IEC 61000-3-2 Class A** standart limitleri tanımlıdır. Analiz edilen cihazın her bir harmoniği limitlerle karşılaştırılarak otomatik **PASS / FAIL** raporu oluşturulur.

---

## 6. Kurulum ve Kullanım

### Gereksinimler
Bu projeyi çalıştırmak için aşağıdaki kütüphanelere sahip Python 3.8 veya daha üstü bir sürüm gereklidir:

```bash
pip install pandas numpy scipy matplotlib opencv-python pillow
