# 1. feladat – A CF-tüdőmikrobiom-hálózat kibővítése és alapparaméterei

> **Helyesbítés a 2. feladat módszertani átdolgozásakor.** Ez a dokumentum az első
> exploratív paraméterkészlet történeti nyilvántartása, **nem validált paraméteradatbázis**.
> A számszerű kölcsönhatások és intervallumaik modellfeltevések; a források sok esetben
> csak asszociációt, más taxont, biofilmet vagy közösségi mechanizmust vizsgálnak.
> A táblákban szereplő „irod.” számok is korábbi feltételezett kiinduló értékek,
> **nem publikációból kiolvasott gLV-együtthatók**. A K-k koncentrációadatok által
> motivált skálapriorok, nem mért eltartóképességek. A növekedési példák egy része
> analógia; az in situ generációs idő nem azonos az intrinsic r-rel.
> A jelenlegi implementációt és e korlátok kezelését a
> [modell_modszertan.md](modell_modszertan.md) írja le. A lenti történeti
> „szakirodalmi/biológiailag ésszerű” jelzők nem jelentenek empirikus validációt.

**Cél:** a tavaszi, 6 fajos Lotka–Volterra modell kibővítése **11 baktériumfajra**, és minden
nem semleges kölcsönhatáshoz, valamint az *r* és *K* paraméterekhez **szakirodalommal alátámasztott
példaérték + finomhangolási intervallum** megadása. A dokumentum a bírálat két fő kritikájára válaszol:
(1) a paraméterek eddig nem voltak forrásokkal indokolva, és „túl szabályosak”, hasraütésszerűek voltak;
(2) az eltartóképesség fizikai dimenziója (absztrakt egység vs. CFU/mL) tisztázatlan volt.

> **Fontos korlát.** A megadott számok **nem egy adott beteg mért paraméterei**, hanem a szakirodalomból
> származtatott, biológiailag ésszerű **példaértékek**, amelyek a megadott intervallumon belül a 3. feladat
> validációs adataival finomhangolhatók. Ahol fajspecifikus mérés nem állt rendelkezésre (jelölve:
> „analógia”), ott rokon faj / ökológiai megfontolás alapján becsült érték szerepel, tág intervallummal.

---

## 1.1. A modell és a mértékegységek

A kiterjesztett, normalizált általánosított Lotka–Volterra (gLV) egyenlet:

$$\frac{dx_i}{dt} = r_i\,x_i\left(1-\frac{x_i}{K_i}+\sum_{j\neq i} a_{ij}\,\frac{x_j}{K_j}\right)$$

| Jelölés | Jelentés | Egység |
|---|---|---|
| *x_i* | az *i* faj abszolút mennyisége | **CFU/mL köpet** |
| *t* | idő | **óra (h)** |
| *r_i* | belső (maximális) növekedési ráta | **1/h** |
| *K_i* | eltartóképesség (hordozókapacitás) | **CFU/mL** |
| *a_ij* | a *j* faj hatása az *i* faj növekedésére | dimenzió nélküli |

**Miért normalizált alak?** A 11 faj *K* értéke több nagyságrendet fog át (10³–10¹⁰ CFU/mL),
ezért a kölcsönhatási tagot a partner saját eltartóképességére normáljuk.
Az *x_j/K_j* nem korlátozott 1-re: facilitáció esetén az egyensúly meghaladhatja K-t.
Az *a_ij* a fajlagos nettó növekedéshez ad *r_i a_ij* hozzájárulást, amikor *x_j=K_j*;
nem általános százalékos populációcsökkenés. Az egyensúlyt az egész hálózat határozza meg.
A normalizálás dimenzionálisan konzisztens, de önmagában nem igazolja az abszolút
CFU-skálát: ehhez mérés vagy explicit kezdeti skálafeltételezés szükséges.

**Konvenció:** `A[i, j]` = a *j*-edik faj (oszlop, ható) hatása az *i*-edik fajra (sor, elszenvedő),
`A[i, i] = −1` (logisztikus önkorlátozás). A hálózati ábrákon a nyíl iránya **j → i**.

---

## 1.2. Fajlista

A tavaszi 6 faj (PA, SA, BCC, SM, Prevotella, Streptococcus) kibővítve a kért 11 fajra. Két funkcionális
csoport: klasszikus **CF-patogének** és **orális/anaerob kommenzálisok** (amelyek a CF-köpetben a
16S-vizsgálatok szerint jelentős részt tesznek ki és befolyásolják a patogének dinamikáját).

| # | Rövidítés | Teljes latin név | Csoport |
|---|---|---|---|
| 1 | PA | *Pseudomonas aeruginosa* | patogén |
| 2 | BCC | *Burkholderia cepacia complex* | patogén |
| 3 | SA | *Staphylococcus aureus* | patogén |
| 4 | HI | *Haemophilus influenzae* | patogén |
| 5 | SM | *Stenotrophomonas maltophilia* | patogén |
| 6 | VD | *Veillonella dispar* | kommenzális (anaerob) |
| 7 | PM | *Prevotella melaninogenica* | kommenzális (anaerob) |
| 8 | SS | *Streptococcus sanguinis* | kommenzális (orális) |
| 9 | FN | *Fusobacterium nucleatum* | kommenzális (anaerob) |
| 10 | RM | *Rothia mucilaginosa* | kommenzális |
| 11 | AX | *Achromobacter xylosoxidans* | patogén |

---

## 1.3. Növekedési ráta (*r*, 1/h)

**Értelmezés.** *r* a belső/maximális (log-fázisú) növekedési ráta: *r* = ln2 / T_kettőződés.
A példaértékeket laboratóriumi (lehetőleg szintetikus CF-köpet médium, SCFM) kettőződési időkből
számítottuk. **Fontos biológiai megállapítás** (Kopf 2016, Kragh 2014): a CF-tüdőben a tényleges
nettó növekedés akár 1–2 nagyságrenddel lassabb, mint in vitro – ez a modellben nem az *r*-ből, hanem a
telítési (1 − *x/K*) tagból és a kölcsönhatásokból adódik. Ezért *r* a felső (in vitro) kapacitást
jelöli, az intervallum alja viszont lefedi a dormáns/in situ állapotot is.

| Rövidítés | Faj | *r* példa [1/h] | *r* intervallum [1/h] | Fő források |
|---|---|---|---|---|
| PA | *Pseudomonas aeruginosa* | 1,3 | 0,05 – 1,66 | Palmer 2007 (SCFM, 32 perc); LaBauve 2012; Kragh 2014; Yang 2011 |
| BCC | *Burkholderia cepacia complex* | 0,5 | 0,02 – 1,0 | analógia (nincs ellenőrzött abszolút érték); Morales 2023 |
| SA | *Staphylococcus aureus* | 0,66 | 0,014 – 1,02 | Yamaguchi 2019 (BHI, medián 63 perc); Kopf 2016 (in situ) |
| HI | *Haemophilus influenzae* | 0,5 | 0,02 – 1,4 | analógia; Do Carmo Silva 2025 (csak relatív r) |
| SM | *Stenotrophomonas maltophilia* | 0,99 | 0,005 – 1,2 | Roscetto 2012 (42 perc); Kolpen 2015 (dormans köpetben) |
| VD | *Veillonella dispar* | 0,5 | 0,01 – 0,8 | Hung 2024 (83 perc, 10 mM laktát; 52 perc nitráttal) |
| PM | *Prevotella melaninogenica* | 0,25 | 0,01 – 0,7 | analógia; Albright 2026 (2-5% O2: fele sebesség) |
| SS | *Streptococcus sanguinis* | 0,42 | 0,02 – 0,6 | Taylor 2025 (TY-glükóz 99 perc; TY-fruktóz 70 perc) |
| FN | *Fusobacterium nucleatum* | 0,2 | 0,005 – 0,3 | Rogers 1991 (kemosztát, µmax 0,20 1/h) |
| RM | *Rothia mucilaginosa* | 0,3 | 0,01 – 0,8 | analógia (nincs ellenőrzött abszolút érték); Gao 2018 |
| AX | *Achromobacter xylosoxidans* | 0,42 | 0,02 – 1,0 | Khademi 2021 (<100 perc = nem lassú izolátum) |

**Számítási példák.** PA: Palmer 2007 SCFM-ben 32 perc kettőződés → *r* = ln2/(32/60) ≈ **1,30 1/h**.
SA: Yamaguchi 2019 BHI-ben medián ~63 perc → ~0,66 1/h; in situ köpetben (Kopf 2016) ~2,1 nap → ~0,014 1/h
(az intervallum alsó vége). SM: Roscetto 2012 42 perc → ~0,99 1/h, de Kolpen 2015 szerint köpetben
gyakorlatilag dormáns (intervallum alja 0,005). FN: Rogers 1991 kemosztát µmax = 0,20 1/h.

---

## 1.4. Eltartóképesség (*K*, CFU/mL)

**Értelmezés és a bírálat kezelése.** A *K* itt **a köpet milliliterére eső maximális
kolóniaképző-egység** – fizikailag jól definiált mennyiség, nem absztrakt egység. Az értékeket
kvantitatív tenyésztéses és qPCR-vizsgálatokból vettük (ahol az adat CFU/g-ban volt, 1 g ≈ 1 mL közelítést
alkalmaztunk). A fajok közti nagyságrendi különbségek **valódiak és irodalmilag alátámasztottak**: pl.
Lamoureux 2021 qPCR-mediánjai szerint Prevotella ~10⁶, Veillonella ~10⁴, Fusobacterium ~10² CFU/mL –
tehát ezek a fajok NEM kezelhetők azonos *K*-val.

| Rövidítés | Faj | *K* példa [CFU/mL] | *K* intervallum [CFU/mL] | Fő források |
|---|---|---|---|---|
| PA | *Pseudomonas aeruginosa* | 5×10^8^ | 10^7^ – 10^10^ | Palmer 2005, 2007; Zemanick 2013 |
| BCC | *Burkholderia cepacia complex* | 10^8^ | 10^6^ – 10^9^ | Carmody 2026; analógia PA-val |
| SA | *Staphylococcus aureus* | 5×10^7^ | 10^5^ – 10^9^ | Worlitzsch 2009; Gangell 2011; Nichols 2023 |
| HI | *Haemophilus influenzae* | 10^7^ | 10^5^ – 10^8^ | Do Carmo Silva 2025; Gangell 2011 |
| SM | *Stenotrophomonas maltophilia* | 10^6^ | 10^4^ – 10^8^ | Kolpen 2015; Terlizzi 2023 |
| VD | *Veillonella dispar* | 10^5^ | 10^3^ – 10^7^ | Lamoureux 2021 (medián 2,7e4) |
| PM | *Prevotella melaninogenica* | 3×10^6^ | 10^4^ – 10^8^ | Lamoureux 2021 (medián 1,3e6); Worlitzsch 2009 |
| SS | *Streptococcus sanguinis* | 10^8^ | 10^6^ – 10^9^ | Sibley 2010 (1e3 - >1e8); Sibley 2011; Worlitzsch 2009 |
| FN | *Fusobacterium nucleatum* | 10^3^ | 10^1^ – 10^5^ | Lamoureux 2021 (medián 2,6e2) |
| RM | *Rothia mucilaginosa* | 10^7^ | 10^5^ – 10^8^ | Rigauts 2022 |
| AX | *Achromobacter xylosoxidans* | 10^7^ | 10^5^ – 10^8^ | Carmody 2026; Menetrey 2021 |

---

## 1.5. Kezdeti sűrűségek (*x₀*, CFU/mL)

A „Teljes terv.docx” példa-mikrobioma (PA 65%, SA 20%, Streptococcus 13%, egyéb 2%) alapján, ~10⁸ CFU/mL
összterhelésre skálázva (felnőtt, krónikus PA-fertőzött beteg). A kisebbségi/kezdetben alig jelen lévő
fajok a detektálási határ közelébe (10³–10⁵) vannak állítva, az anaerobok a Lamoureux 2021 mediánok körül.

| Rövidítés | Faj | *x₀* [CFU/mL] | Megjegyzés |
|---|---|---|---|
| PA | *Pseudomonas aeruginosa* | 6,5×10^7^ | 65 % |
| BCC | *Burkholderia cepacia complex* | 10^5^ | kisebbségi |
| SA | *Staphylococcus aureus* | 2×10^7^ | 20 % |
| HI | *Haemophilus influenzae* | 10^5^ | kisebbségi |
| SM | *Stenotrophomonas maltophilia* | 10^5^ | kisebbségi |
| VD | *Veillonella dispar* | 3×10^4^ | Lamoureux 2021 medián |
| PM | *Prevotella melaninogenica* | 10^6^ | Lamoureux 2021 medián |
| SS | *Streptococcus sanguinis* | 1,3×10^7^ | 13 % |
| FN | *Fusobacterium nucleatum* | 3×10^2^ | Lamoureux 2021 medián |
| RM | *Rothia mucilaginosa* | 5×10^5^ | kisebbségi |
| AX | *Achromobacter xylosoxidans* | 10^5^ | kisebbségi |

---

## 1.6. Kölcsönhatási mátrix (*a_ij*)

**76 nem semleges, irányított kapcsolat** (30 gátló, 46 támogató) – a tavaszi 6×6-os mátrix ~10
kapcsolatához képest. Minden kapcsolat legalább egy elsődleges forrással alátámasztott; az erősségi
besorolás (gyenge/közepes/erős) a bizonyíték típusából és mértékéből adódik:

- **erős** (|*a*| ≈ 0,6–0,9): közvetlen ölés / >1 log CFU-változás (pl. PA→SA HQNO-val, PA→BCC);
- **közepes** (|*a*| ≈ 0,3–0,6): 50–90% hatás vagy több egybehangzó tanulmány;
- **gyenge** (|*a*| ≈ 0,05–0,3): <50% hatás vagy csak együttes-előfordulási (ko-okkurrencia) korreláció.

**A „példaérték” és az „irodalmi besorolás” viszonya.** A példaértékeket az irodalmi besorolásból indítva,
**a szakirodalmi intervallumon belül** finomítottuk (lásd 1.7.), hogy antibiotikum nélkül mind a 11 faj
stabil együttélést mutasson. Ahol a finomított érték eltér az irodalmi besorolástól, az utóbbit a táblázat
*(irod. …)* jelöléssel feltünteti; a kódban (paraméter-cella) ugyanez szerepel megjegyzésként. Az intervallum
minden esetben a szakirodalmi erősségi tartomány.

> **Előjel-konfliktusok.** Néhány kapcsolatnál a laboratóriumi és a klinikai ko-okkurrencia bizonyíték
> ellentétes (pl. PA↔SS: in vitro PA serkentheti az orális streptococcusokat, de a CF-hálózatokban
> negatívan korrelálnak; SA↔HI: Margolis 2010 +19% vs. Pettigrew 2008 negatív korreláció). Ezeknél az
> intervallum **átnyúlik a nullán**, és a példaérték a konzervatívabb/klinikailag relevánsabb irányt követi.

| # | Ható *j* | Elszenvedő *i* | Jel | *a_ij* példa | Intervallum | Fő források |
|---|---|---|---|---|---|---|
| 1 | PA *Pseudomonas aeruginosa* | SA *Staphylococcus aureus* | 🟥 gátló | -0,55 *(irod. -0,75)* | -0,90 … -0,30 | Hoffman 2006; Filkins 2015; Limoli 2017 (mukoid PA: koegzisztencia) |
| 2 | SA *Staphylococcus aureus* | PA *Pseudomonas aeruginosa* | 🟩 támogató | +0,20 | +0,05 … +0,40 | Mashburn 2005; Choi 2024 |
| 3 | PA *Pseudomonas aeruginosa* | BCC *Burkholderia cepacia complex* | 🟥 gátló | -0,65 *(irod. -0,70)* | -0,90 … -0,50 | Bragonzi 2012; Tomlin 2001; Bakkal 2010 |
| 4 | BCC *Burkholderia cepacia complex* | PA *Pseudomonas aeruginosa* | 🟩 támogató | +0,15 | -0,30 … +0,40 | Bragonzi 2012 (biofilm) vs. Bakkal 2010 (bakteriocin) |
| 5 | PA *Pseudomonas aeruginosa* | HI *Haemophilus influenzae* | 🟥 gátló | -0,15 | -0,30 … -0,05 | Quinn 2016; Lyczak 2002 |
| 6 | HI *Haemophilus influenzae* | PA *Pseudomonas aeruginosa* | 🟩 támogató | +0,10 | +0,00 … +0,25 | Duan 2003 |
| 7 | PA *Pseudomonas aeruginosa* | SM *Stenotrophomonas maltophilia* | 🟥 gátló | -0,65 *(irod. -0,70)* | -0,90 … -0,50 | Pompilio 2015; Menetrey 2020 |
| 8 | SM *Stenotrophomonas maltophilia* | PA *Pseudomonas aeruginosa* | 🟩 támogató | +0,20 | +0,05 … +0,40 | Ryan 2008; Pompilio 2015; Magalhaes 2016 |
| 9 | PA *Pseudomonas aeruginosa* | AX *Achromobacter xylosoxidans* | 🟥 gátló | -0,20 *(irod. -0,10)* | -0,30 … +0,00 | Menetrey 2020; Sandri 2021 |
| 10 | AX *Achromobacter xylosoxidans* | PA *Pseudomonas aeruginosa* | 🟥 gátló | -0,55 | -0,75 … -0,35 | Menetrey 2020; Sandri 2021 |
| 11 | PA *Pseudomonas aeruginosa* | VD *Veillonella dispar* | 🟥 gátló | -0,25 *(irod. -0,15)* | -0,30 … -0,05 | Quinn 2016 |
| 12 | VD *Veillonella dispar* | PA *Pseudomonas aeruginosa* | 🟩 támogató | +0,10 | +0,00 … +0,20 | Flynn 2016; Mirkovic 2015 |
| 13 | PA *Pseudomonas aeruginosa* | PM *Prevotella melaninogenica* | 🟥 gátló | -0,40 *(irod. -0,30)* | -0,50 … -0,10 | Quinn 2016 |
| 14 | PM *Prevotella melaninogenica* | PA *Pseudomonas aeruginosa* | 🟩 támogató | +0,25 | +0,10 … +0,45 | Flynn 2016; Moore 2026; Sherrard 2016 |
| 15 | PA *Pseudomonas aeruginosa* | SS *Streptococcus sanguinis* | 🟥 gátló | -0,10 *(irod. +0,15)* | -0,10 … +0,40 | Whiley 2015 (in vitro +) vs. Quinn 2016 |
| 16 | SS *Streptococcus sanguinis* | PA *Pseudomonas aeruginosa* | 🟥 gátló | -0,50 | -0,90 … -0,20 | Scoffield & Wu 2015, 2016 |
| 17 | PA *Pseudomonas aeruginosa* | FN *Fusobacterium nucleatum* | 🟥 gátló | -0,10 | -0,25 … -0,05 | Quinn 2016 |
| 18 | FN *Fusobacterium nucleatum* | PA *Pseudomonas aeruginosa* | 🟩 támogató | +0,30 | +0,10 … +0,50 | Flynn 2016 |
| 19 | PA *Pseudomonas aeruginosa* | RM *Rothia mucilaginosa* | 🟥 gátló | -0,40 *(irod. -0,30)* | -0,50 … -0,10 | Vandeplassche 2020; Quinn 2016; Layeghifard 2019 |
| 20 | RM *Rothia mucilaginosa* | PA *Pseudomonas aeruginosa* | 🟩 támogató | +0,20 | +0,05 … +0,40 | Gao 2018 |
| 21 | SA *Staphylococcus aureus* | BCC *Burkholderia cepacia complex* | 🟩 támogató | +0,15 | +0,00 … +0,40 | Bernier 2016 |
| 22 | SM *Stenotrophomonas maltophilia* | BCC *Burkholderia cepacia complex* | 🟩 támogató | +0,10 | +0,00 … +0,30 | Bernier 2016 |
| 23 | SA *Staphylococcus aureus* | HI *Haemophilus influenzae* | 🟩 támogató | +0,15 | -0,30 … +0,30 | Margolis 2010 (+19%) vs. Pettigrew 2008 |
| 24 | HI *Haemophilus influenzae* | SA *Staphylococcus aureus* | 🟥 gátló | -0,05 *(irod. -0,10)* | -0,30 … +0,00 | Pettigrew 2008 |
| 25 | SA *Staphylococcus aureus* | SM *Stenotrophomonas maltophilia* | 🟥 gátló | -0,05 | -0,15 … +0,00 | Menetrey 2020 |
| 26 | AX *Achromobacter xylosoxidans* | SM *Stenotrophomonas maltophilia* | 🟥 gátló | -0,25 *(irod. -0,30)* | -0,50 … -0,10 | Menetrey 2020, 2021 |
| 27 | SS *Streptococcus sanguinis* | SA *Staphylococcus aureus* | 🟥 gátló | -0,30 *(irod. -0,35)* | -0,50 … -0,20 | Oogai 2016 |
| 28 | FN *Fusobacterium nucleatum* | SA *Staphylococcus aureus* | 🟩 támogató | +0,35 *(irod. +0,30)* | +0,10 … +0,50 | Lima 2019 |
| 29 | SA *Staphylococcus aureus* | FN *Fusobacterium nucleatum* | 🟩 támogató | +0,30 | +0,10 … +0,50 | Lima 2019 |
| 30 | VD *Veillonella dispar* | SA *Staphylococcus aureus* | 🟥 gátló | -0,05 *(irod. -0,10)* | -0,30 … +0,00 | Quinn 2016 |
| 31 | SA *Staphylococcus aureus* | VD *Veillonella dispar* | 🟥 gátló | -0,15 *(irod. -0,10)* | -0,30 … +0,00 | Quinn 2016 |
| 32 | PM *Prevotella melaninogenica* | SA *Staphylococcus aureus* | 🟥 gátló | -0,15 | -0,30 … +0,00 | Quinn 2016; Muhlebach 2018 |
| 33 | SA *Staphylococcus aureus* | PM *Prevotella melaninogenica* | 🟥 gátló | -0,15 | -0,30 … +0,00 | Quinn 2016; Muhlebach 2018 |
| 34 | RM *Rothia mucilaginosa* | SA *Staphylococcus aureus* | 🟥 gátló | -0,05 *(irod. -0,10)* | -0,30 … +0,00 | Muhlebach 2018 |
| 35 | SS *Streptococcus sanguinis* | HI *Haemophilus influenzae* | 🟥 gátló | -0,15 | -0,30 … -0,05 | Pericone 2000 (analógia); Kreth 2005 |
| 36 | SS *Streptococcus sanguinis* | VD *Veillonella dispar* | 🟩 támogató | +0,30 *(irod. +0,50)* | +0,30 … +0,80 | Egland 2004; Periasamy 2010; Chalmers 2008 |
| 37 | VD *Veillonella dispar* | SS *Streptococcus sanguinis* | 🟩 támogató | +0,20 *(irod. +0,40)* | +0,20 … +0,60 | Egland 2004; Periasamy 2010 |
| 38 | FN *Fusobacterium nucleatum* | SS *Streptococcus sanguinis* | 🟩 támogató | +0,10 *(irod. +0,30)* | +0,10 … +0,60 | Kaplan 2009; Periasamy 2010 |
| 39 | SS *Streptococcus sanguinis* | FN *Fusobacterium nucleatum* | 🟥 gátló | -0,45 *(irod. -0,20)* | -0,50 … +0,20 | Zhou 2017 |
| 40 | VD *Veillonella dispar* | FN *Fusobacterium nucleatum* | 🟩 támogató | +0,35 *(irod. +0,50)* | +0,30 … +0,80 | Zhou 2017; Periasamy 2010 |
| 41 | FN *Fusobacterium nucleatum* | VD *Veillonella dispar* | 🟩 támogató | +0,25 *(irod. +0,40)* | +0,20 … +0,60 | Periasamy 2010 |
| 42 | PM *Prevotella melaninogenica* | FN *Fusobacterium nucleatum* | 🟩 támogató | +0,25 *(irod. +0,30)* | +0,10 … +0,50 | Okuda 2012 |
| 43 | FN *Fusobacterium nucleatum* | PM *Prevotella melaninogenica* | 🟩 támogató | +0,20 *(irod. +0,30)* | +0,10 … +0,50 | Okuda 2012 |
| 44 | VD *Veillonella dispar* | PM *Prevotella melaninogenica* | 🟩 támogató | +0,10 *(irod. +0,15)* | +0,05 … +0,30 | Flynn 2016; Quinn 2016; Segal 2013 |
| 45 | PM *Prevotella melaninogenica* | VD *Veillonella dispar* | 🟩 támogató | +0,10 *(irod. +0,15)* | +0,05 … +0,30 | Flynn 2016; Quinn 2016; Segal 2013 |
| 46 | PM *Prevotella melaninogenica* | SS *Streptococcus sanguinis* | 🟩 támogató | +0,15 *(irod. +0,20)* | +0,10 … +0,35 | Flynn 2016; Quinn 2016 |
| 47 | SS *Streptococcus sanguinis* | PM *Prevotella melaninogenica* | 🟩 támogató | +0,15 *(irod. +0,20)* | +0,10 … +0,35 | Flynn 2016; Quinn 2016 |
| 48 | SS *Streptococcus sanguinis* | RM *Rothia mucilaginosa* | 🟩 támogató | +0,20 | +0,10 … +0,30 | Quinn 2016; Layeghifard 2019 |
| 49 | RM *Rothia mucilaginosa* | SS *Streptococcus sanguinis* | 🟩 támogató | +0,10 *(irod. +0,20)* | +0,10 … +0,30 | Quinn 2016; Layeghifard 2019 |
| 50 | FN *Fusobacterium nucleatum* | RM *Rothia mucilaginosa* | 🟩 támogató | +0,05 *(irod. +0,10)* | +0,05 … +0,30 | Quinn 2016 |
| 51 | RM *Rothia mucilaginosa* | FN *Fusobacterium nucleatum* | 🟩 támogató | +0,05 *(irod. +0,10)* | +0,05 … +0,30 | Quinn 2016 |
| 52 | VD *Veillonella dispar* | RM *Rothia mucilaginosa* | 🟩 támogató | +0,05 *(irod. +0,10)* | +0,05 … +0,30 | Quinn 2016 |
| 53 | RM *Rothia mucilaginosa* | VD *Veillonella dispar* | 🟩 támogató | +0,05 *(irod. +0,10)* | +0,05 … +0,30 | Quinn 2016 |
| 54 | PM *Prevotella melaninogenica* | RM *Rothia mucilaginosa* | 🟩 támogató | +0,10 | +0,05 … +0,30 | Quinn 2016 |
| 55 | RM *Rothia mucilaginosa* | PM *Prevotella melaninogenica* | 🟩 támogató | +0,05 *(irod. +0,10)* | +0,05 … +0,30 | Quinn 2016 |
| 56 | VD *Veillonella dispar* | BCC *Burkholderia cepacia complex* | 🟩 támogató | +0,10 | +0,00 … +0,25 | Flynn 2016 |
| 57 | PM *Prevotella melaninogenica* | BCC *Burkholderia cepacia complex* | 🟩 támogató | +0,10 | +0,00 … +0,25 | Flynn 2016 |
| 58 | SS *Streptococcus sanguinis* | BCC *Burkholderia cepacia complex* | 🟩 támogató | +0,25 *(irod. +0,20)* | +0,05 … +0,35 | Flynn 2016 |
| 59 | FN *Fusobacterium nucleatum* | BCC *Burkholderia cepacia complex* | 🟩 támogató | +0,15 | +0,05 … +0,30 | Flynn 2016 |
| 60 | BCC *Burkholderia cepacia complex* | VD *Veillonella dispar* | 🟥 gátló | -0,20 *(irod. -0,10)* | -0,25 … +0,00 | Muhlebach 2018; Coburn 2015 |
| 61 | BCC *Burkholderia cepacia complex* | PM *Prevotella melaninogenica* | 🟥 gátló | -0,40 *(irod. -0,30)* | -0,50 … -0,10 | Cuthbertson 2020; Muhlebach 2018 |
| 62 | BCC *Burkholderia cepacia complex* | SS *Streptococcus sanguinis* | 🟥 gátló | -0,35 *(irod. -0,20)* | -0,35 … -0,05 | Coburn 2015; Flight 2015 |
| 63 | BCC *Burkholderia cepacia complex* | FN *Fusobacterium nucleatum* | 🟥 gátló | -0,15 *(irod. -0,10)* | -0,25 … +0,00 | Muhlebach 2018 |
| 64 | VD *Veillonella dispar* | SM *Stenotrophomonas maltophilia* | 🟩 támogató | +0,10 | +0,00 … +0,25 | Flynn 2016 |
| 65 | PM *Prevotella melaninogenica* | SM *Stenotrophomonas maltophilia* | 🟩 támogató | +0,10 | +0,00 … +0,25 | Flynn 2016 |
| 66 | SS *Streptococcus sanguinis* | SM *Stenotrophomonas maltophilia* | 🟩 támogató | +0,30 *(irod. +0,25)* | +0,10 … +0,40 | Flynn 2016; Muhlebach 2018 |
| 67 | FN *Fusobacterium nucleatum* | SM *Stenotrophomonas maltophilia* | 🟩 támogató | +0,15 *(irod. +0,10)* | -0,05 … +0,30 | Flynn 2016 vs. Muhlebach 2018 |
| 68 | SM *Stenotrophomonas maltophilia* | VD *Veillonella dispar* | 🟥 gátló | -0,15 *(irod. -0,10)* | -0,25 … +0,00 | Muhlebach 2018 |
| 69 | SM *Stenotrophomonas maltophilia* | PM *Prevotella melaninogenica* | 🟥 gátló | -0,15 *(irod. -0,10)* | -0,30 … +0,05 | Cuthbertson 2020; Muhlebach 2018 |
| 70 | SM *Stenotrophomonas maltophilia* | SS *Streptococcus sanguinis* | 🟩 támogató | +0,05 *(irod. +0,15)* | +0,00 … +0,30 | Muhlebach 2018 |
| 71 | SM *Stenotrophomonas maltophilia* | FN *Fusobacterium nucleatum* | 🟥 gátló | -0,20 *(irod. -0,15)* | -0,30 … +0,00 | Muhlebach 2018 |
| 72 | SM *Stenotrophomonas maltophilia* | RM *Rothia mucilaginosa* | 🟩 támogató | +0,10 *(irod. +0,15)* | +0,00 … +0,30 | Muhlebach 2018 |
| 73 | VD *Veillonella dispar* | AX *Achromobacter xylosoxidans* | 🟩 támogató | +0,05 *(irod. +0,10)* | +0,00 … +0,25 | Flynn 2016 |
| 74 | PM *Prevotella melaninogenica* | AX *Achromobacter xylosoxidans* | 🟩 támogató | +0,05 *(irod. +0,10)* | +0,00 … +0,25 | Flynn 2016 |
| 75 | SS *Streptococcus sanguinis* | AX *Achromobacter xylosoxidans* | 🟩 támogató | +0,10 *(irod. +0,20)* | +0,05 … +0,35 | Flynn 2016 |
| 76 | FN *Fusobacterium nucleatum* | AX *Achromobacter xylosoxidans* | 🟩 támogató | +0,10 *(irod. +0,15)* | +0,05 … +0,30 | Flynn 2016 |

---

## 1.7. Numerikus együttélési hangolás – nem adatkalibráció

Az első verzióban r és K rögzítése mellett a 76 feltételezett kölcsönhatásból 49
értékét numerikusan módosítottuk az előzetesen választott intervallumokon belül.
A cél pozitív, lokálisan stabil együttélési egyensúly volt, nem mért adatok illesztése.
A 0,2–1 közötti normalizált egyensúly csak büntetőfüggvényes cél volt, nem teljesült
minden komponensre kemény korlátként. A ~66%-os PA-arány ebből a választásból ered;
nem független validációs eredmény.

A korábbi notebook csak az egyensúly ellenőrzését tartalmazta, nem az optimalizálót.
A 2. feladat nem optimalizálja újra ezeket az értékeket: pontos megőrzésüket regressziós
teszt ellenőrzi, és kikapcsolható, hipotetikus hálózatként használja őket.
Pozitív egyensúly lokális stabilitása nem garantál globális stabilitást vagy
biológiai hitelességet, különösen nem kezelés alatt és más paraméterekkel.

---

## 1.8. Ábrák

Az 1. feladatból megőrzött ábrák (a 2. feladat dinamikai kimenetei az `eredmenyek`
könyvtárban készülnek, nem ezeket az archív ábrákat írják felül):

| Fájl | Tartalom |
|---|---|
| `mikrobiom_halozat.png` | Teljes kölcsönhatási hálózat (11 csomópont, 76 nyíl; piros=gátló, zöld=támogató, vastagság ∝ |*a*|) |
| `mikrobiom_halozat_eros_kapcsolatok.png` | Letisztult hálózat, csak |*a*| ≥ 0,25 |
| `kolcsonhatasi_matrix.png` | Az *a_ij* mátrix hőtérképe |
| `r_K_parameterek.png` | *r* és *K* példaérték + finomhangolási intervallum (log skálán) |
| `lotka_volterra_dinamika_log.png`* | Szimulált log₁₀(CFU/mL) dinamika (a notebookban jelenik meg) |

(*a dinamikaábrák a notebook kimeneteként láthatók.)

---

## 1.9. Kimenetek (a feladat szerint)

1. **Táblázatos, forrásokkal hivatkozott paraméterek** (*r*, *K*, kölcsönhatási mátrix) – ez a dokumentum.
2. **Mikrobiomhálózat-ábra** – `mikrobiom_halozat.png` (+ letisztult változat).
3. **Kibővített `lotka_volterra_alap.ipynb`** – 11 faj, normalizált CFU/mL-es gLV, szakirodalmi paraméterek,
   hálózati és dinamikai ábrák, egyensúly-ellenőrzés.

---

## 1.10. Módszertani megjegyzések és korlátok

- **A kölcsönhatások nagy része ko-okkurrencia-alapú** (Quinn 2016, Muhlebach 2018, Layeghifard 2019,
  Cuthbertson 2020), ami korreláció, nem ok-okozat. A legmegbízhatóbb, mechanizmussal is bizonyított
  élek: PA→SA (HQNO), SS→PA (H₂O₂/nitrit), FN↔SA/SS/VD (koaggregáció), SS↔VD (laktát-keresztfeedelés),
  PM→PA (β-laktamáz-védelem).
- **Az AX és RM kapcsolatai a leggyengébben alátámasztottak**; Cuthbertson 2020 kifejezetten kivételként
  kezeli az AX-et. Ezek tág intervallummal, gyenge példaértékkel szerepelnek.
- **BCC/SM/AX ↔ kommenzális élek** jórészt egyetlen keresztfeedelési kísérletre (Flynn 2016, konzorcium-
  szinten) és egyetlen faktoranalízisre (Muhlebach 2018) támaszkodnak – a bizalom alacsony-közepes.
- Az intervallumok célja épp az, hogy a 3. feladat validációs adatai (TIP/S, AZLI, azithromycin
  mikrobiom-idősorok) alapján a paraméterek **adatvezérelten finomhangolhatók** legyenek.

---

## 1.11. Irodalomjegyzék

A hivatkozások a fenti táblázatokban szerzőnévvel és évvel szerepelnek. Minden forrást a kutatás során
PubMed / PMC / Europe PMC teljes szövegből vagy absztraktból ellenőriztünk (a DOI/PMCID valós).

### Növekedési ráta és eltartóképesség (*r*, *K*)

1. Palmer KL, Aye LM, Whiteley M. (2007) *Nutritional cues control Pseudomonas aeruginosa multicellular behavior in cystic fibrosis sputum.* J Bacteriol 189(22):8079–8087. doi:10.1128/JB.01138-07. PMC2168676.
2. Palmer KL, Mashburn LM, Singh PK, Whiteley M. (2005) *Cystic fibrosis sputum supports growth and cues key aspects of Pseudomonas aeruginosa physiology.* J Bacteriol 187(15):5267–5277. PMC1196007.
3. LaBauve AE, Wargo MJ. (2012) *Growth and laboratory maintenance of Pseudomonas aeruginosa.* Curr Protoc Microbiol. PMC4296558.
4. Turner KH, Wessel AK, Palmer GC, Murray JL, Whiteley M. (2015) *Essential genome of Pseudomonas aeruginosa in cystic fibrosis sputum.* PNAS 112(13):4110–4115. PMC4386324.
5. Kragh KN, et al. (2014) *Polymorphonuclear leukocytes restrict growth of Pseudomonas aeruginosa in the lungs of cystic fibrosis patients.* Infect Immun 82(11):4477–4486. PMC4249348.
6. Kopf SH, et al. (2016) *Trace incorporation of heavy water reveals slow and heterogeneous pathogen growth rates in cystic fibrosis sputum.* PNAS 113(2):E110–E116. PMC4720290.
7. Yang L, et al. (2011) *Evolutionary dynamics of bacteria in a human host environment.* PNAS 108(18):7481–7486. PMC3088582.
8. Yamaguchi T, et al. (2019) *Association between cell growth and vancomycin resistance in clinical community-associated MRSA.* Infect Drug Resist 12:2379–2390. PMC6682180.
9. Roscetto E, et al. (2012) *Functional characterization of the RNA chaperone Hfq in Stenotrophomonas maltophilia.* J Bacteriol 194(21):5864–5874. PMC3486118.
10. Kolpen M, et al. (2015) *Denitrification by cystic fibrosis pathogens – Stenotrophomonas maltophilia is dormant in sputum.* Int J Med Microbiol 305(1):1–10. PMID 25441256.
11. Morales LD, Av-Gay Y, Murphy MEP. (2023) *Acidic pH modulates Burkholderia cenocepacia antimicrobial susceptibility in the CF nutritional environment.* Microbiol Spectr 11(6):e0273123. PMC10714822.
12. Do Carmo e Silva P, Hill D, Harrison F. (2025) *Optimizing synthetic cystic fibrosis sputum media for growth of non-typeable Haemophilus influenzae.* Access Microbiol. PMC12181625.
13. Hung J-H, Zhang S-M, Huang S-L. (2024) *Nitrate promotes the growth and the production of SCFAs and tryptophan from Veillonella dispar.* Appl Environ Microbiol 90(8):e01148-24. PMC11337843.
14. Albright C, et al. (2026) *Aerotolerant capacity of the lung symbiont Prevotella melaninogenica.* J Bacteriol 208(8):e00142-26. PMC13492135.
15. Taylor ZA, Pham DN, Zeng L. (2025) *Systematic analysis of the glucose-PTS in Streptococcus sanguinis…* Appl Environ Microbiol 91(1):e01935-24. PMC11784104.
16. Rogers AH, et al. (1991) *Aspects of the growth and metabolism of Fusobacterium nucleatum ATCC 10953 in continuous culture.* Oral Microbiol Immunol 6(4):250–255. PMID 1812468.
17. Khademi SMH, et al. (2021) *[Achromobacter xylosoxidans CF clinical isolates characterization].* mSystems. PMC8269239.
18. Gao B, et al. (2018) *Tracking Polymicrobial Metabolism in CF Airways: Pseudomonas aeruginosa … influenced by Rothia mucilaginosa-Derived Metabolites.* mSphere 3(2):e00151-18. PMC5917424.
19. Zemanick ET, et al. (2013) *Inflammation and airway microbiota during cystic fibrosis pulmonary exacerbations.* PLoS One 8(4):e62917. PMC3639911.
20. Tunney MM, et al. (2008) *Detection of anaerobic bacteria in high numbers in sputum from patients with cystic fibrosis.* Am J Respir Crit Care Med 177(9):995–1001. PMID 18263800.
21. Worlitzsch D, et al. (2009) *Antibiotic-resistant obligate anaerobes during exacerbations of cystic fibrosis patients.* Clin Microbiol Infect 15(5):454–460. PMID 19196263.
22. Lamoureux C, et al. (2021) *An observational study of anaerobic bacteria in cystic fibrosis lung…* Sci Rep 11:6845. PMC7994387.
23. Gangell C, et al. (2011) *Inflammatory responses to individual microorganisms in the lungs of children with CF.* Clin Infect Dis 53(5):425–432. PMID 21844026.
24. Sibley CD, et al. (2010) *McKay agar enables routine quantification of the 'Streptococcus milleri' group in CF patients.* J Med Microbiol 59:534–540. PMID 20093379.
25. Sibley CD, et al. (2011) *Culture enriched molecular profiling of the cystic fibrosis airway microbiome.* PLoS One 6(7):e22702. PMC3145661.
26. Rigauts C, et al. (2022) *Rothia mucilaginosa is an anti-inflammatory bacterium in the respiratory tract of patients with chronic lung disease.* Eur Respir J 59(5):2101293. PMC9068977.
27. Carmody LA, et al. (2026) *Density fluctuations of Burkholderia and Achromobacter in cystic fibrosis sputum.* J Cyst Fibros 25(4):651–654. PMID 41832128.
28. Nichols DP, et al. (2023) *Pharmacologic improvement of CFTR function rapidly decreases sputum pathogen density…* J Clin Invest 133(10):e167957. PMC10178839.
29. Menetrey Q, et al. (2021) *Achromobacter xylosoxidans and Stenotrophomonas maltophilia: Emerging Pathogens…* Genes 12(5):610. PMC8142972.
30. Terlizzi V, et al. (2023) *Stenotrophomonas maltophilia in people with Cystic Fibrosis: a systematic review…* Eur J Clin Microbiol Infect Dis 42(11):1285–1296. PMC10587323.

### Kölcsönhatások

31. Hoffman LR, et al. (2006) *Selection for Staphylococcus aureus small-colony variants due to growth in the presence of Pseudomonas aeruginosa.* PNAS 103(52):19890–19895. PMC1750898.
32. Filkins LM, et al. (2015) *Coculture of Staphylococcus aureus with Pseudomonas aeruginosa drives S. aureus towards fermentative metabolism…* J Bacteriol 197(14):2252–2264. PMC4524177.
33. Limoli DH, et al. (2017) *Pseudomonas aeruginosa Alginate Overproduction Promotes Coexistence with Staphylococcus aureus…* mBio 8(2):e00186-17. PMC5362032.
34. Mashburn LM, et al. (2005) *Staphylococcus aureus serves as an iron source for Pseudomonas aeruginosa during in vivo coculture.* J Bacteriol 187(2):554–566. PMC543556.
35. Choi SY, et al. (2024) *Autolysis of Pseudomonas aeruginosa Quorum-Sensing Mutant Is Suppressed by Staphylococcus aureus through Iron-Dependent Metabolism.* J Microbiol Biotechnol 34(4):795–803. PMC11091660.
36. Bragonzi A, et al. (2012) *Modelling co-infection of the CF lung by Pseudomonas aeruginosa and Burkholderia cenocepacia…* PLoS One 7(12):e52330. PMC3528780.
37. Tomlin KL, Coll OP, Ceri H. (2001) *Interspecies biofilms of Pseudomonas aeruginosa and Burkholderia cepacia.* Can J Microbiol 47(10):949–954. PMID 11718549.
38. Bakkal S, et al. (2010) *Role of bacteriocins in mediating interactions of bacterial isolates taken from cystic fibrosis patients.* Microbiology 156:2058–2067. PMC3068677.
39. Lyczak JB, Cannon CL, Pier GB. (2002) *Lung infections associated with cystic fibrosis.* Clin Microbiol Rev 15(2):194–222. PMC118069.
40. Duan K, et al. (2003) *Modulation of Pseudomonas aeruginosa gene expression by host microflora through interspecies communication.* Mol Microbiol 50(5):1477–1491. PMID 14651632.
41. Pompilio A, et al. (2015) *Cooperative pathogenicity in cystic fibrosis: Stenotrophomonas maltophilia modulates Pseudomonas aeruginosa virulence in mixed biofilm.* Front Microbiol 6:951. PMC4584994.
42. Ryan RP, et al. (2008) *Interspecies signalling via the Stenotrophomonas maltophilia diffusible signal factor influences biofilm formation and polymyxin tolerance in Pseudomonas aeruginosa.* Mol Microbiol 68(1):75–86. PMID 18312265.
43. Magalhães AP, Lopes SP, Pereira MO. (2016) *Insights into Cystic Fibrosis Polymicrobial Consortia…* Front Microbiol 7:2146. PMC5233685.
44. Menetrey Q, et al. (2020) *High Occurrence of Bacterial Competition Among Clinically Documented Opportunistic Pathogens Including Achromobacter xylosoxidans in CF.* Front Microbiol 11:558160. PMC7513574.
45. Sandri A, et al. (2021) *Adaptive Interactions of Achromobacter spp. with Pseudomonas aeruginosa in CF Chronic Lung Co-Infection.* Pathogens 10(8):978. PMC8400197.
46. Flynn JM, Niccum D, Dunitz JM, Hunter RC. (2016) *Evidence and Role for Bacterial Mucin Degradation in Cystic Fibrosis Airway Disease.* PLoS Pathog 12(8):e1005846. PMC4993466.
47. Mirković B, et al. (2015) *The Role of Short-Chain Fatty Acids, Produced by Anaerobic Bacteria, in the Cystic Fibrosis Airway.* Am J Respir Crit Care Med 192(11):1314–1324. PMC4731701.
48. Moore PJ, Kent LA, Hunter RC. (2026) *Anaerobic microbiota promote pathogen association with the airway epithelium.* J Med Microbiol 75(3):002149. PMC13016583.
49. Sherrard LJ, et al. (2016) *Production of extended-spectrum β-lactamases and the potential indirect pathogenic role of Prevotella isolates from the CF respiratory microbiota.* Int J Antimicrob Agents 47(2):140–145. PMC4746055.
50. Whiley RA, et al. (2015) *Environment and colonisation sequence are key parameters driving cooperation and competition between Pseudomonas aeruginosa CF strains and oral commensal streptococci.* PLoS One 10(2):e0115513. PMC4339374.
51. Scoffield JA, Wu H. (2015) *Oral streptococci and nitrite-mediated interference of Pseudomonas aeruginosa.* Infect Immun 83(1):101–107. PMC4288860.
52. Scoffield JA, Wu H. (2016) *Nitrite reductase is critical for Pseudomonas aeruginosa survival during co-infection with Streptococcus parasanguinis.* Microbiology 162(2):376–383. PMC4766596.
53. Vandeplassche E, et al. (2020) *Antibiotic susceptibility of cystic fibrosis lung microbiome members in a multispecies biofilm.* Biofilm 2:100031. PMC7798459.
54. Layeghifard M, et al. (2019) *Microbiome networks and change-point analysis reveal key community changes associated with cystic fibrosis pulmonary exacerbations.* npj Biofilms Microbiomes 5:4. PMC6341074.
55. Bernier SP, et al. (2016) *Cyanide Toxicity to Burkholderia cenocepacia Is Modulated by Polymicrobial Communities and Environmental Factors.* Front Microbiol 7:725. PMC4870242.
56. Margolis E, Yates A, Levin BR. (2010) *The ecology of nasal colonization of Streptococcus pneumoniae, Haemophilus influenzae and Staphylococcus aureus…* BMC Microbiol 10:59. PMC2844402.
57. Pettigrew MM, et al. (2008) *Microbial interactions during upper respiratory tract infections.* Emerg Infect Dis 14(10):1584–1591. PMC2609881.
58. Oogai Y, Kawada-Matsuo M, Komatsuzawa H. (2016) *Staphylococcus aureus SrrAB Affects Susceptibility to Hydrogen Peroxide and Co-Existence with Streptococcus sanguinis.* PLoS One 11(7):e0159768. PMC4956065.
59. Lima BP, et al. (2019) *The Oral Bacterium Fusobacterium nucleatum Binds Staphylococcus aureus and Alters Expression of the Staphylococcal Accessory Regulator sarA.* Microb Ecol 78(2):336–347. PMC11796435.
60. Quinn RA, et al. (2016) *Ecological networking of cystic fibrosis lung infections.* npj Biofilms Microbiomes 2:16004 (2:4). doi:10.1038/s41522-016-0002-1. PMC5460249. (a kiindulási példaforrás)
61. Muhlebach MS, et al. (2018) *Anaerobic bacteria cultured from cystic fibrosis airways correlate to milder disease: a multisite study.* Eur Respir J 52(1):1800242. PMC6376871.
62. Pericone CD, et al. (2000) *Inhibitory and bactericidal effects of hydrogen peroxide production by Streptococcus pneumoniae on other inhabitants of the upper respiratory tract.* Infect Immun 68(7):3990–3997. PMC101678. (fajanalógia S. pneumoniae → viridans streptococcus)
63. Kreth J, et al. (2005) *Competition and coexistence between Streptococcus mutans and Streptococcus sanguinis in the dental biofilm.* J Bacteriol 187(21):7193–7203. PMC1272965.
64. Egland PG, Palmer RJ Jr, Kolenbrander PE. (2004) *Interspecies communication in Streptococcus gordonii–Veillonella atypica biofilms…* PNAS 101(48):16917–16922. PMC534724.
65. Periasamy S, Kolenbrander PE. (2010) *Central role of the early colonizer Veillonella sp. in establishing multispecies biofilm communities…* J Bacteriol 192(12):2965–2972. PMC2901697.
66. Chalmers NI, et al. (2008) *Characterization of a Streptococcus sp.–Veillonella sp. community micromanipulated from dental plaque.* J Bacteriol 190(24):8145–8154. PMC2593232.
67. Kaplan CW, et al. (2009) *The Fusobacterium nucleatum outer membrane protein RadD is an arginine-inhibitable adhesin…* Mol Microbiol 71(1):35–47. PMC2741168.
68. Zhou P, et al. (2017) *Veillonella catalase protects the growth of Fusobacterium nucleatum in microaerophilic and Streptococcus gordonii-resident environments.* Appl Environ Microbiol 83(19):e01079-17. PMC5601340.
69. Okuda T, et al. (2012) *Synergy in biofilm formation between Fusobacterium nucleatum and Prevotella species.* Anaerobe 18(1):110–116. PMID 21945879.
70. Segal LN, et al. (2013) *Enrichment of lung microbiome with supraglottic taxa is associated with increased pulmonary inflammation.* Microbiome 1:19. PMC3971609.
71. Cuthbertson L, et al. (2020) *Lung function and microbiota diversity in cystic fibrosis.* Microbiome 8(1):45. PMC7114784.
72. Coburn B, et al. (2015) *Lung microbiota across age and disease stage in cystic fibrosis.* Sci Rep 5:10241. PMC4431465.
73. Flight WG, et al. (2015) *Rapid Detection of Emerging Pathogens and Loss of Microbial Diversity Associated with Severe Lung Disease in Cystic Fibrosis.* J Clin Microbiol 53(7):2022–2029. PMC4473198.

> **A kutatás részletes, idézetekkel és mennyiségi adatokkal ellátott jegyzőkönyvei** (háttéranyag,
> nem a leadandó dokumentum része) a session `files/` mappájában: `research_rK.md`, `research_rK_addendum.md`,
> `research_growth_rates.md`, `research_PA.md`, `research_pathogens.md`, `research_commensals.md`,
> `research_gramneg_commensal.md`.
