# 2. feladat – A CF-mikrobiommodell matematikai átdolgozása

## 1. Státusz és értelmezési határ

**Ez tesztelhető kutatási modellkeret, nem klinikailag validált betegmodell és nem
adagolási tanácsadó.** Az itt megadott példafutás nem igazolja a kezelés hatásosságát,
biztonságát vagy a rezisztencia várható klinikai időpontját. A numerikus verifikáció,
a paraméterkalibráció és a független biológiai validáció három külön feladat.

Az átdolgozás alapjai:

- a [bírálattal ellátott dolgozat](TDK_2026_Schleier_Anna_biralat.pdf) tényleges PDF-megjegyzései,
  különösen a 14–16., 19., 24., 27., 29–32. és 36. oldalon;
- a [módszertani fejlesztési terv](TDK_dolgozat_fejlesztesi_es_atdolgozasi_terv.docx),
  valamint az abban megadott Spalding-, Nashebi- és Landersdorfer-cikkek [1–4];
- a [CFU–relatív abundancia dokumentum](CFU_és_relatív_abundancia_kapcsolata.pdf).
  Ez AI-válasz, nem elsődleges szakirodalom: a kompozicionalitás kérdését [9–11] is alátámasztja;
- a felhasználó által megerősített döntések: **napi összdózis mg-ban; állandó abszolút
  Others-háttér; demonstrációs immunprofil; gamma-bemenet + elimináció, külön PD**.

### Mit javít a megvalósítás?

| Korábbi probléma | Új megoldás |
|---|---|
| Antibiotikumok önkényes r/K-szorzói, akár negatív r vagy K | Változatlan pozitív K; nemnegatív növekedésgátlás és külön koncentrációfüggő ölés |
| Gauss-hatásablak és külön exponenciális lecsengés | Kauzális gamma/Erlang bejutás, helyi elimináció; nincs második önkényes hatásablak |
| Dózis dimenziótlan, nincs helyi koncentrációskála | mg/dózis → becsült helyi mg/L, explicit átváltási paraméterrel |
| Korlátlan, lineárisan változó R-szint; kiválasztott gyógyszerek hibás indexelése | Gyógyszerenként külön rezisztenciabitek, tényleges szubpopulációk, fix gyógyszersorrend |
| Általános, nem azonosítható HGT | Nincs HGT; szelekció és születéshez kötött mutáció |
| Kemény „kihalási” küszöb az ODE-ban | Nincs ilyen küszöb; a detektálási határ csak megfigyelési réteg |
| Csak kiválasztott taxonok újranormalizálása | Explicit Others-háttér; a megjelenített teljes összetétel 100% |
| Paraméterek és egyenletek egy notebookban | Fix modell, szerkeszthető paraméterek, külön futtató |
| Stabil szimuláció biológiai validációként értelmezve | Forrásadat, átvitt hipotézis és numerikus teszt egyértelmű szétválasztása |

## 2. Fájlok és futtatás

1. **[mikrobiom_modell.py](mikrobiom_modell.py)**: egyenletek, adatszerkezetek,
   inputellenőrzés, dózisesemények, integrálás, megfigyelési függvények.
   Paraméterillesztéskor ezt nem kell módosítani.
2. **[modell_parameterek.py](modell_parameterek.py)**: 11 taxon, r/K és interakciós
   táblák, PK/PD, immunprofil, kezdeti mikrobiom, napi kezelési sorok.
3. **[szimulacio_futtatas.py](szimulacio_futtatas.py)**: Python API és parancssori
   futtatás, ábrák és adatexport.
4. **[lotka_volterra_alap.ipynb](lotka_volterra_alap.ipynb)**: vékony interaktív
   belépési pont, ugyanazokat a modulokat használja, nem tart fenn második ODE-t.
5. **[test_mikrobiom.py](test_mikrobiom.py)**: matematikai és regressziós tesztek.

A korábbi notebook megőrzött példánya:
[lotka_volterra_alap_safety_copy.ipynb](lotka_volterra_alap_safety_copy.ipynb).
A hálózati számok változatlan átvitelét ez ellenőrzi; a régi antibiotikum-modell
hibáit nem kell megőrizni. A korábbi PNG-k archív Task-1 eredmények.

A projekt könyvtárában, kiválasztott Python-környezetben:

```powershell
python .\szimulacio_futtatas.py
python .\szimulacio_futtatas.py --snapshots 0 2.5 7 10 14
python .\szimulacio_futtatas.py --no-antibiotics --output .\eredmenyek_kontroll
python .\szimulacio_futtatas.py --regimen-csv .\kezeles.csv
python -m unittest test_mikrobiom -v
```

Ezen a gépen a `python` rendszeralias helyett szükség lehet erre:

```powershell
& 'C:\Users\z005awpa\AppData\Local\Programs\Python\Python314\python.exe' .\szimulacio_futtatas.py
```

Futási függőségek: NumPy, SciPy, Matplotlib; a notebookhoz IPython/Jupyter,
az interaktív hálózati ábrához NetworkX. Nincs futás közbeni automatikus
csomagtelepítés. A környezet verzióit a futás JSON-ja elmenti.

### Egyszerű szerkesztési sorrend

- `SELECTED_SPECIES`: tetszőleges, nem üres részhalmaz a 11 rövidítésből.
- `INITIAL_MODE` és a megfelelő kezdeti értékek: a kulcsok pontosan egyezzenek a
  kiválasztott fajokkal; az ismeretlen taxon vagy hiányzó PD-sor hiba.
- `PATIENT`: nem, életkor, BMI, tüdőkategória, társbetegségek.
- `REGIMEN`: naponta, antibiotikumonként egy sor. Üres tuple: nincs kezelés.
- `DURATION_DAYS`, `SNAPSHOT_DAYS`, `OUTPUT_STEP_H`, `MAX_STEP_H`.
- Illesztéskor `R_TABLE`, `K_TABLE`, `INTERACTIONS`, `PK_PARAMETERS`,
  `PD_PARAMETERS`. A PD-számokat az alap `PD_ROWS` táblában vagy a létrejövő
  `PD_PARAMETERS` egyedi rekordjainak cseréjével lehet módosítani.

A `PD` és `PK` rekordok immutábilisak; az egyedi változtatás mintája:

```python
from dataclasses import replace
import modell_parameterek as p

p.PD_PARAMETERS["PA"]["tobramycin"] = replace(
    p.PD_PARAMETERS["PA"]["tobramycin"], ec50_mg_l=80.0
)
model = p.build_model(("PA", "SA"))
```

A később létrehozott `Model` kapja az új értékeket; egy korábban összerakott
modell újraépítése szükséges. Az r/K/A példákat a kód a hozzájuk megadott
intervallummal is ellenőrzi; a PK/PD tuningintervallumok útmutatók, a fizikai
és matematikai tartományokat a rekordok külön validálják.

## 3. Idő, dózis és mértékegységek

| Mennyiség | Egység és konvenció |
|---|---|
| Belső idő t | óra |
| X-tengely | t/24, nap |
| Kezelési nap | 1-indexelt; 1. nap = [0,24), 8. nap = [168,192) h |
| Szimuláció vége | 14 nap = 336 h; ezen az időponton már nincs új dózis |
| Napi adag | napi **összes** mg, nem mg/beadás |
| Beadásonkénti adag | napi_mg / beadásszám |
| Tranzitkompartment Q | mg |
| Helyi koncentráció C | mg/L = µg/mL |
| Baktérium N | CFU/mL; a példainicializálásnál becsült CFU/mL |
| r, ölési ráta, immun-clearance | h⁻¹ |
| K | CFU/mL; környezetfüggő modellparaméter |

Ha a napi beadásórák nincsenek megadva, a kód egyenletesen osztja el őket:
`0, 24/n, ..., 24(n−1)/n`. Ez kényelmi **szimulációs konvenció**, nem
klinikai adagolási előírás. Az órák külön is megadhatók.

```python
from mikrobiom_modell import DailyDose

kezeles = (
    DailyDose("tobramycin", day=1, daily_mg=600, administrations=2, hours=(8, 20)),
    DailyDose("aztreonam", day=1, daily_mg=750, administrations=3),
    DailyDose("azithromycin", day=8, daily_mg=300, administrations=2),
)
```

CSV-alternatíva:

```csv
drug,day,daily_mg,administrations,hours
aztreonam,1,750,3,0;8;16
aztreonam,2,0,0,
azithromycin,8,300,2,8;20
```

A hiányzó napok kezelésmentesek. Ugyanarra a gyógyszerre és napra két sor,
negatív dózis, ismételt beadásóra, horizonton kívüli beadás explicit hibát okoz.
0 mg-os sorhoz 0 beadás szükséges. A CLI CSV csak a kezelést cseréli; a beteg,
fajok és időhorizont továbbra is a paramétermodulból jönnek.

**A felhasználói példafutás:** aztreonam 750 mg/nap, 3×250 mg, 1–14. nap;
azithromycin 300 mg/nap, 2×150 mg, 8–14. nap; tobramycin kikapcsolva.
Az előbbi aztreonamadag eltér a Cayston forrásban szereplő 75 mg naponta
háromszor sémától [6]. Nem javítottuk át önkényesen, de nem is ajánljuk.
A megemelt dózisra alkalmazott lineáris expozícióskálázás nincs validálva.

## 4. Farmakokinetika: gamma-bemenet és helyi kiürülés

### 4.1. Egyenletek

Gyógyszerenként m egész számú tranzitkompartmentet használunk; a példában m=2.
Ez a gamma-eloszlás egész alakparaméterű, **Erlang** esete, nem minden lehetséges
nem egész gamma-alak megvalósítása. Legyen:

\[
k_{\rm tr}=m/\tau,\qquad k_{\rm el}=\ln(2)/T_{1/2,\rm local}.
\]

A t_d dózisidőben csak az első tranzitállapot ugrik:

\[
Q_1(t_d^+)=Q_1(t_d^-)+D_d,\qquad D_d\ [\mathrm{mg}].
\]

Dózisok között:

\[
\dot Q_1=-k_{\rm tr}Q_1,\qquad
\dot Q_\ell=k_{\rm tr}(Q_{\ell-1}-Q_\ell),\quad \ell=2,\ldots,m,
\]
\[
\dot C=\beta\,k_{\rm tr}Q_m-k_{\rm el}C.
\]

Itt \(\beta\) egysége L⁻¹, a beadott mg és a modell által reprezentált
helyi mg/L közötti **effektív dózis–expozíció tényező**.

Egy dózis helyi bejutási sebessége:

\[
k_{\rm tr}Q_m(t)=D_d f_\Gamma(t-t_d;m,k_{\rm tr}),\qquad
f_\Gamma(u)=\frac{k_{\rm tr}^{m}u^{m-1}e^{-k_{\rm tr}u}}{\Gamma(m)}
\mathbf{1}_{u\ge0}.
\]

A gamma-bemenet integrálja D_d. A lokális koncentráció:

\[
C(t)=\beta\sum_{t_d\le t}D_d
\int_0^{t-t_d} f_\Gamma(u)e^{-k_{\rm el}(t-t_d-u)}\,du.
\]

**Minden korábbi dózis maradéka automatikusan összeadódik**: nincs állapotnullázás
napváltáskor vagy új beadáskor. Az egyetlen dózis teljes AUC-ja:
\(\mathrm{AUC}_{0,\infty}=\beta D/k_{\rm el}\), egysége mg·h/L.
A kezdeti Q és C nulla: előzetes krónikus gyógyszerexpozíciót a példafutás nem
feltételez. Az első dózis előtt nincs hatás, és nem történik közvetlen C-ugrás.

### 4.2. Mit nem lehet az orális/inhalált dózisból egyedül megmondani?

A mg nem koncentráció. A \(\beta=F_{\rm local}/V_{\rm effective}\) felbontás csak akkor
lenne fizikailag azonosítható, ha a helyi hozzáférhetőség és effektív térfogat is
ismert lenne. Itt csak összevont \(\beta\)-t adunk meg, nem kitalált „tüdőtérfogatot”.

- Inhalációnál a készülék, aeroszol, belégzéstechnika és lerakódás befolyásolja a
  koncentrációt. A névleges dózis nem teljes egészében jut a mintázott régióba.
- Orális azithromycinnél az abszorpció, szöveti és sejten belüli felhalmozódás,
  visszaáramlás is lényeges [7,8]. A modell ezt effektív tranzit- és
  eliminációs paraméterekbe sűríti; **nem teljes test-PK**.
- Az ábra a feltételezett helyi koncentrációt mutatja, nem a teljes szervezetben
  maradt mg-ot vagy a szérumkoncentrációt. Utóbbihoz külön kompartmentek és adatok kellenek.
- Köpet teljes gyógyszertartalom ≠ szabad, biofilmben ható koncentráció.
  A PD-paramétereket a választott helyi mérési skálához kell illeszteni.
- A µg/g köpetadat mg/L-é alakításához explicit sűrűségfeltételezés kell:
  \(C_{\rm mg/L}=C_{\rm \mu g/g}\rho_{\rm g/mL}\). A példában \(\rho=1\);
  ez nem a CFU/g és CFU/mL minden esetben pontos azonossága.

### 4.3. Mért forrásadat és implementált prior külön táblázatban

| Szer | Ellenőrzött forrásadat | Mit nem bizonyít? |
|---|---|---|
| Inhalált tobramycin | Geller 2002: 300 mg után átlagos köpetcsúcs 1237 µg/g; 1 h-nál szérum 0,95 µg/mL; becsült szisztémás biohasznosulás 11,7%; vizsgálati sémában nincs felhalmozódás [5] | Nem ad univerzális helyi eliminációt vagy minden eszközre dózisarányos expozíciót |
| Inhalált aztreonam | EMA: 75 mg után 10 perccel átlagos köpet 726 µg/g; 1 h-nál plazma ~0,59 µg/mL; **szérum** eliminációs felezési idő ~2,1 h [6] | A 2,1 h nem köpetfelezési idő; a 726 egy mintavételi érték, nem teljes PK-görbe |
| Orális azithromycin | Wilms 2006, n=8 CF felnőtt, hosszú távú 500 mg/nap: köpet 12–53 mg/L; 10 nappal leállítás után 4–27 mg/L; plazma felezési idő 102±20 h, PMNN 289±166 h [7] | Ezek nem egyetlen dózis helyi csúcsai, nem közvetlen paraméterei a mi 14 napos kezdő sémánknak |
| Orális azithromycin | SmPC: formulafüggő biohasznosulás, plazma Tmax 2–3 h; terminális plazmafelezési idő ~68 h [8] | A plazma Tmax/felezési idő nem azonos a köpetével |

| Paraméter | Tobramycin példa [tuning] | Aztreonam példa [tuning] | Azithromycin példa [tuning] |
|---|---|---|---|
| m, gamma-alak | 2 | 2 | 2 |
| τ, átlagos tranzitidő [h] | 0,10 [0,02–0,5] | 0,10 [0,02–0,5] | 2,0 [0,5–6] |
| T½, helyi [h] | 1,0 [0,25–4] | 1,0 [0,25–4] | 72 [24–168] |
| β [L⁻¹] | 5,2 [0,5–20] | 12,2 [1–40] | 0,010 [0,001–0,05] |

**Mindegyik implementált helyi PK-szám modellprior, nem betegadatból illesztett
becslés.** Az inhalált β értékei a 10 perces köpetkoncentrációk nagyságrendjét
közelítik a feltételezett tranzit/felezési idő mellett; ez egy ponthoz igazítás,
nem kinetikai validáció. Az azithromycin-prior hosszabb helyi jelenlétet reprezentál,
de a 72 h nem a [7]-ből kiolvasott köpetfelezési idő. A tuningintervallumok
exploratív keresési tartományok, nem konfidenciaintervallumok.

## 5. Ökológia, antibiotikumhatás és rezisztencia együtt

### 5.1. Miért nyolc fenotípus fajonként?

A három gyógyszerhez három bit tartozik, állandó sorrendben:
**tobramycin, aztreonam, azithromycin**.
A fenotípusok: `000, 001, 010, 011, 100, 101, 110, 111`.
A 0 a modellben érzékenyebb, az 1 az adott szerre kevésbé érzékeny fenotípus.
Ez nem klinikai breakpoint szerinti érzékenységi diagnózis.

Legyen \(X_{ib}\ge0\) az i taxon b fenotípusának CFU/mL-je és
\(N_i=\sum_bX_{ib}\). Három külön S/R pár ugyanannak a fajnak a sejtjeit
többször számolná; egyetlen „mindenre rezisztens R” pedig indokolatlan
keresztrezisztenciát feltételezne. A nyolc állapot mindkettőt elkerüli.
Egyes fenotípusok pontosan nullák maradhatnak; nincs mesterséges pozitív vetés.

A kezdeti gyógyszerspecifikus frakciókból a joint eloszlás:
\[
X_{ib}(0)=N_i(0)\prod_{j=1}^3 f_{ij}^{b_j}(1-f_{ij})^{1-b_j}.
\]
A kódban ezt elágazással számoljuk, így f=0 és f=1 is megengedett.
**Ez kezdeti függetlenségi feltevés**, nem mért multirezisztens klóneloszlás.
Ha ilyen adatok rendelkezésre állnak, az inicializálást kell arra lecserélni.

### 5.2. gLV előjel és nemnegatív születés/veszteség felbontás

Megmarad az első feladat konvenciója:
\[
\dot N_i=r_iN_i\left(1-\frac{N_i}{K_i}+
\sum_{j\ne i}a_{ij}\frac{N_j}{K_j}\right).
\]
`A[i,j]` a j faj i fajra gyakorolt hatása; diagonális −1.
Az ábra nyila j→i, pozitív támogató, negatív gátló.

Definiáljuk \(a^+_{ij}=\max(a_{ij},0)\), \(a^-_{ij}=\max(-a_{ij},0)\), j≠i:
\[
B_i(N)=r_i\left(1+\sum_{j\ne i}a^+_{ij}N_j/K_j\right),
\qquad
L_i(N)=r_i\left(N_i/K_i+\sum_{j\ne i}a^-_{ij}N_j/K_j\right)+I.
\]
Így B és L nemnegatív. A bakteriostatikus hatás csak a pozitív B születési
tagot szorozza, **nem a negatív nettó gLV-deriváltat**. Ez elkerüli, hogy
egy növekedésgátló a negatív nettó növekedés „gátlásával” paradox módon
megmentse a populációt. A születés/veszteség felbontás ugyanazt a nettó
gyógyszermentes gLV-t adja, de **új mechanisztikus feltevés**: a nettó
görbe önmagában nem azonosítja a sejtosztódási és pusztulási fluxusokat.

### 5.3. Hill/Emax PD, külön ölés és növekedésgátlás

\[
E_{ibj}(C_j)=\frac{C_j^{h_{ij}}}
{(\mathrm{EC50}_{ij}F_{ij}^{b_j})^{h_{ij}}+C_j^{h_{ij}}}.
\]
F a rezisztens EC50-szorzó, F≥1. EC50 nem azonos automatikusan a MIC-kel.
Az azonos EC50 a kill/inhibition modulban tudatos minimumparaméterezés;
a jelenlegi nem nulla sorok csak az egyik mechanizmust kapcsolják be.

\[
g_{ib}(C)=\prod_j(1-I_{\max,ij}E_{ibj}),\qquad
k_{ib}(C)=\sum_j k_{\max,ij}E_{ibj},
\qquad
w_{ib}=\prod_j(1-c_{ij})^{b_j}.
\]

0≤Imax≤1, kmax≥0, 0≤c≤1, ezért a gyógyszer nem állít elő negatív
növekedési faktort vagy negatív ölési rátát. Tobramycin és aztreonam:
külön baktericid veszteség; azithromycin: a példában bakteriostatikus
növekedésgátlás az explicit kijelölt taxonokon [1,4,6,8,12].

**Kombinációs alapfeltevés:** az ölési hazardok összege, a növekedési
megmaradó faktorok szorzata. Nem implementáltunk forrás nélküli szinergia-
vagy antagonizmusszámokat. A Nashebi-cikk koncepcionális S/R mintát ad,
nem a mi három szerünk illesztett kombinációs paramétereit. Több szer
egyidejű futtatása működik, de nagyobb összhatás nem bizonyít szinergiát.
A kész keret jelenleg **nem kalibrál explicit gyógyszer–gyógyszer
interakciós paramétert**; ez kombinációs time-kill adatokkal külön bővítés.

### 5.4. Mutáció és a teljes ODE

A j szerrel szemben még érzékeny fenotípusból az egy bitben eltérő rezisztensbe:
\[
J_{ibj}=\mu_{ij}\,w_{ib}g_{ib}B_i\,X_{ib},\qquad b_j=0.
\]
\(\mu\) dimenziótlan osztódásonkénti valószínűségi prior; nincs
koncentrációküszöbös „rezisztencianövelés”. Fajonként \(\sum_j\mu_{ij}\le1\).
Egy eseményben egy rezisztenciabit változik, további bitek későbbi
átmenetekkel szerezhetők; nincs közvetlen több-bites mutáció.

\[
\boxed{
\dot X_{ib}=
\left[w_{ib}g_{ib}B_i-L_i-k_{ib}\right]X_{ib}
+\sum_{j:b_j=1}J_{i,b-e_j,j}
-\sum_{j:b_j=0}J_{ibj}
}
\]

A mutációs tagok fajon belül összegezve kiesnek: sejtszámot
nem teremtenek és nem semmisítenek meg. A rezisztenciaarány:
\[
R_{ij}(t)=\frac{\sum_{b:b_j=1}X_{ib}(t)}{N_i(t)}.
\]
N_i=0 esetén a frakció nem definiált, a kód NaN-t, a CSV üres mezőt
használ, nem hamis 0%-ot. A modell érvényes állapotainál egyébként [0,1].

Nincs visszamutáció, adaptív toleranciaállapot, HGT, immunsejtállapot
vagy biofilmkompartment. Ezek hiánya **nem annak állítása, hogy biológiailag
nem léteznek**, hanem az azonosíthatatlan bővítések elkerülése.

### 5.5. Gyorsan „kialakul” vagy szelekcióval felszaporodik a rezisztencia?

Mindkettő lehetséges, és a mechanizmusok elkülönítendők:

1. **Előzetes rezisztens kisebbség szelekciója:** már f(0)>0; az érzékenyek
   gyorsabban fogynak, ezért az R-frakció mutáció nélkül is nőhet.
2. **De novo genetikai változás:** a fenti J-fluxus után a leszármazottak
   szelekcióval felszaporodhatnak. A mutációs gyakoriság nem azonos az
   osztódásonkénti mutációs rátával.
3. **Reverzibilis adaptív rezisztencia/tolerancia:** például aminoglikozid-
   adaptáció [4]; ezt a jelenlegi genetikai bitmodell nem írja le.

Landersdorfer 2018-ban a hiper-mutábilis PAOΔmutS HFIM-monoterápiáiban
a kevésbé érzékeny baktériumok a 3. napra közel teljesen felválthatták az
eredeti populációt; PAO1 mellett is volt gyors kezdeti ölés utáni
visszanövekedés [3]. Ezek konkrét törzsek, expozíciók és in vitro rendszerek,
**nem univerzális „3 nap alatt rezisztenssé válik a beteg” szabályok**.
A meropenem–tobramycin eredmények nem az aztreonam–tobramycin kombináció
kvantitatív validációi.

## 6. PD-paraméterek és bizonyítéki státuszuk

**Az alábbi Emax/EC50 számok mind nem illesztett demonstrációs priorszámok.**
A források a mechanizmust vagy az aktivitásspektrum egy részét támasztják alá,
nem a számokat. BCC komplex, nem egyetlen egységes érzékenységű faj.
Az orális taxonokra azithromycinhatás átvitele különösen gyenge feltevés.

Jelölések a táblában: `K(kmax; EC50)` = ölés [h⁻¹; mg/L],
`S(Imax; EC50)` = növekedésgátlás [1; mg/L], `0` = explicit
nullhipotézis / nem modellezett közvetlen hatás, **nem bizonyított semlegesség**.

| Taxon | Tobramycin inhalált | Aztreonam inhalált | Azithromycin orális | Forrás és átvitel korlátja |
|---|---|---|---|---|
| P. aeruginosa | K(2,0;100) | K(1,8;100) | 0 | Tobramycin és aztreonam mechanizmus [3–6]; azithromycin CF-előnyét nem azonosítjuk közvetlen PA-öléssel [7,8] |
| B. cepacia complex | 0 | 0 | 0 | Tobramycin insusceptible lista [4]; más szerekre adat nélküli nullhipotézis, nem komplexszintű általános állítás |
| S. aureus | K(1,0;100) | 0 | S(0,85;1) | Tobramycin válaszképes listán [4]; aztreonam Gram-negatív spektrum [6]; makrolidaktivitás rezisztenciafüggő [8] |
| H. influenzae | K(1,0;100) | K(1,0;100) | S(0,85;1) | Tobramycin [4]; aztreonam Gram-negatív spektrumából exploratív átvitel [6]; makrolidspektrum [8] |
| S. maltophilia | 0 | 0 | 0 | Tobramycin insusceptible lista [4]; más szereknél nullhipotézis |
| V. dispar | 0 | 0 | 0 | Nincs beillesztett faj-/törzsspecifikus time-kill adat |
| P. melaninogenica | 0 | 0 | 0 | Nincs beillesztett faj-/törzsspecifikus time-kill adat |
| S. sanguinis | 0 | 0 | S(0,80;1) | Makrolidmechanizmusból [8] átvitt, nem fajspecifikusan validált demonstráció |
| F. nucleatum | 0 | 0 | 0 | Nincs beillesztett faj-/törzsspecifikus time-kill adat |
| R. mucilaginosa | 0 | 0 | S(0,80;1) | Hipotetikus makrolidválasz; nincs ebből fajspecifikus Emax/EC50 adat |
| A. xylosoxidans | 0 | 0 | 0 | Tobramycin insusceptible lista régi Alcaligenes névvel [4]; egyéb közvetlen hatások nincsenek modellezve |

A nullsorok explicit szerepelnek a kódban, ezért egy kihagyott faj–szer
rekord nem válik észrevétlenül „hatástalanná”.
Minden sornál Hill=2, F=16. Nem nulla hatású soroknál:
f(0)=10⁻⁴, μ=10⁻⁸/osztódás, c=0,05. Nullsoroknál f=μ=c=0.
Ezek a rezisztenciapriorok **nem** [3] paramétereinek átmásolásai.

| Paraméter | Példa | Tuningtartomány | Bizonyítéki státusz |
|---|---|---|---|
| kmax [h⁻¹] | fenti sorok | 0–5 | time-killből illesztendő; függvényalak [1,12] |
| Imax | 0–0,85 | 0–1 | demonstráció; növekedési görbéből illesztendő |
| EC50 [mg/L] | 1 vagy 100 | 0,01–1000 | skála- és törzsfüggő prior; MIC-ből nem következik közvetlenül |
| Hill h | 2 | 1–4 | prior |
| Rezisztens EC50-szorzó F | 16 | 1–1000 | prior; rezisztens time-kill szükséges |
| Kezdeti f | 10⁻⁴ | 0–0,1 | prior; szelektív tenyésztési adat kell |
| μ/osztódás | 10⁻⁸ | 0–10⁻⁵ | prior; mutációs vizsgálat kell |
| Fitneszköltség c | 0,05 | 0–0,3 | prior; gyógyszermentes S/R-növekedés kell |

Kombinációs kölcsönhatás, mutáció és költség nélkül, azonos fenotípus-PD esetén
visszakapjuk a megfelelő egyszerűsített modellhatárt. A költség a születési
fluxust csökkenti; a változatlan K ellenére a rezisztens monokultúra
egyensúlya ettől alacsonyabb lehet, mert a veszteségi fluxus nem csökken vele.

## 7. r, K és a 76 kölcsönhatás újraértelmezése

A [Task-1 paramétertáblák](mikrobiom_parameterek.md) r/K/A értékeit és
intervallumait számszerű változtatás nélkül vittük át. Ez összehasonlíthatóság,
**nem a korábbi állítások megerősítése**.

- Ko-előfordulásból nem következik irányított kauzális gLV-együttható.
- Biofilmtömeg, koaggregáció, virulencia és nettó CFU-növekedés külön végpont.
- Konzorciumhatás nem feltétlenül bontható minden tagra külön páronkénti élként.
- K nem a köpetben mért medián koncentráció vagy a tenyésztési kimutatási küszöb.
- In vitro kettőződési idő, in situ turnover és intrinsic r nem felcserélhető.
- A korábbi stabil együttélésre történő numerikus hangolás nem adatillesztés.
- Az intervallumon belüli tetszőleges kombináció **nem garantáltan stabil**.

Az `INTERACTION_SCALE=0` kontroll kikapcsolja a fajok közti éleket, az
önkorlátozást nem. Az 1-es példa a korábbi hipotetikus hálózat.
Validálás előtt indokolt ritkább, valódi együtt-tenyésztési adatokkal
azonosítható hálózatot is összevetni vele.

Gyógyszer, mutáció, fitneszköltség és immunitás nélkül a pozitív egyensúly:
\[
y^*=-A^{-1}\mathbf1,\quad N^*=\mathrm{diag}(K)y^*,\quad
J_y=\mathrm{diag}(r_i y_i^*)A.
\]
Csak y*>0 esetén értelmezhető all-positive egyensúlyként; lokálisan stabil,
ha minden sajátérték valós része negatív. A notebook ezt a **tiszta ökológiai
határesetet** ellenőrzi, nem a kezelt fenotípusrendszer egyensúlyát.
Kanonikus gLV alakban a dimenziós keresztegyüttható
\(\alpha_{ij}=r_i a_{ij}/K_j\), egysége mL/(CFU·h).

## 8. Kezdeti CFU és relatív abundancia

### 8.1. Egyetlen abszolút horgony

Adott kezdeti p_i százalékok és az a horgony N_a koncentrációja:
\[
N_{\rm total}(0)=\frac{N_a(0)}{p_a/100},\qquad
N_i(0)=N_a(0)\frac{p_i}{p_a},\qquad
O=N_{\rm total}(0)\left(1-\sum_i p_i/100\right).
\]

A példában N_PA=6,5·10⁷ CFU/mL, p_PA=65%, p_SA=20%, p_SS=13%:
N_total=10⁸, N_SA=2·10⁷, N_SS=1,3·10⁷, O=2·10⁶ CFU/mL.
A 65/20/13/2 összetétel így pontos, nem közelítően összeállított.
Az abszolút horgony **példa**, nem mérés; az export is becsültként jelöli.

Ez a képlet matematikailag helyes, ha a p_i arány ugyanannak az abszolút
mennyiségnek az aránya. 16S-read arány és életképes CFU-arány általában
nem azonos: eltérő génkópiaszám, DNS-extrakció, PCR, halott sejtek,
tenyészthetőség, telepképző aggregátumok torzíthatnak [9–11].
Ideálisan a referencia ugyanabból a betegből, mintából és időpontból származik.
Az így becsült értékek **nem független abszolút validációs mérések**.
Egyetlen t=0 horgony nem határozza meg a későbbi abszolút terhelést.

Alternatívák: `log10` módban \(N_i=10^{z_i}\); `cfu` módban közvetlen N_i.
Pontos nulla csak a közvetlen CFU bemenetben vagy 0%-ként adható meg,
mert log10(0) nem véges. A log10 bemenetben az `Others` külön abszolút adat.

### 8.2. Others a szimuláció alatt

A felhasználó választása szerint \(\dot O=0\):
\[
p_i(t)=100\frac{N_i(t)}{O+\sum_jN_j(t)},\qquad
p_{\rm Others}(t)=100\frac{O}{O+\sum_jN_j(t)}
=100-\sum_i p_i(t).
\]

Az Others aránya változik, abszolút mennyisége nem. Nincs rá antibiotikum,
immun-clearance vagy ökológiai kölcsönhatás. Ez egyszerű, explicit háttérfeltevés,
nem a kihagyott taxonok tényleges modellje. Részhalmaz váltásakor a felhasználónak
újra meg kell adnia az összetételt: a kód nem teszi automatikusan az összes
elhagyott fajt az Others kategóriába.

### 8.3. Megfigyelés és későbbi illesztés

A kimeneti logmegfigyelés:
\[
z_i(t)=\log_{10}\max(N_i(t),L_{\rm det}).
\]
A példa L_det=1 CFU/mL, **szerkeszthető demonstrációs határ**, nem a mintavétel
kimutatási határának megállapítása. Ez a csonkolás csak a megjelenítést érinti;
a nyers CFU és fenotípusok megmaradnak. A Δlog panel a cenzorált z(t)−z(0).
Nem jelent bizonyított sterilizációt, ha egy görbe eléri a padlót.

A relatív 16S-idősorokat elsősorban a saját mérési skálájukon kell illeszteni.
Az alap megfigyelési modell p_i∝N_i, de nem kalibrált 16S-átalakítás.
Később fajonkénti torzítási súlyokkal \(p_i^{seq}\propto b_iN_i\)
és megfelelő kompozicionális hibamodellel bővíthető; ez **nincs most illesztve**.
Kizárólag relatív adatokból az abszolút K-skála általában nem azonosítható.

## 9. Pácienskategóriákból számított immunkonstans

Nem, életkor és BMI ellenőrzött, exportált metaadat; a dinamika és a PK
számításába nem kerülnek. A kategóriák:

| Tüdőkategória | FEV1 | Demonstrációs w |
|---|---|---|
| mild | ≥70% | 1,0 |
| moderate | ≥40 és <70% | 0,7 |
| severe | <40% | 0,4 |

Kategóriát kell megadni, nem a FEV1-ből számolunk rejtett többletparamétert.
Az immun-clearance:
\[
I=I_0w_{\rm lung}\exp\left(-\sum_{d\in\mathrm{társbetegségek}}v_d\right),
\quad I_0=0{,}01\ \mathrm{h}^{-1}.
\]

| Társbetegség | Konfigurációs kulcs | v példa |
|---|---|---|
| Pancreatic insufficiency | pancreatic_insufficiency | 0,05 |
| CF-related diabetes | cf_related_diabetes | 0,15 |
| Impaired glucose tolerance | impaired_glucose_tolerance | 0,08 |
| Liver disease | liver_disease | 0,05 |
| Bone disease | bone_disease | 0,02 |
| Sinus disease | sinus_disease | 0,02 |
| Recurrent DIOS | recurrent_dios | 0,03 |

**A képlet, I0 és minden súly a felhasználó által kért demonstrációs modellfeltevés.
Nincs hozzá validált CF-specifikus átváltási formula.** FEV1 nem immunmérés,
és a CF gyulladás fokozottsága nem azonos a hatékony bakteriális clearance-szel.
A súlyok csökkenő clearance-t modelleznek, nem bizonyított betegség-okozati hatást.

Feltáró tartomány: I0=0–0,05 h⁻¹, w=0–1, v=0–0,5; ezek is szabadon választott
forgatókönyv-tartományok. Kontrollhoz I0=0; kategóriák semlegesítéséhez minden
w=1 és minden v=0. Diabetes és impaired glucose tolerance egyszerre nem
jelölhető, mert ebben az egyszerű kategóriamodellben kettős számolást okozna.
A példaprofil moderate + pancreatic insufficiency:
I=0,01·0,7·exp(−0,05)≈0,0066586 h⁻¹.

## 10. Numerikus eljárás és verifikáció

- `solve_ivp(method="LSODA")`: az inhalált bejutás gyors, az azithromycin
  helyi jelenléte lassú, ezért merevségre is alkalmas solver.
- A horizont minden dózisidőnél szegmensekre bomlik; ugrás kizárólag Q1-ben.
  Nem lép át a solver észrevétlenül egy dózisdiszkontinuitást.
- Példa: `max_step_h=0.5`, `output_step_h=0.5`, rtol=10⁻⁷, atol=10⁻⁹.
  Az adaptív belső lépés lehet 0,1 h-nál kisebb, ha a gyors tranzit igényli.
  A 0,1–1 h nem biztonságos kötelező alsó lépéshatár.
- Minden kért pillanat és dózisidő bekerül a kimeneti rácsba.
  A gyors gamma-felfutás láthatóságához a dózis utáni 0,5τ, τ, 2τ, 4τ
  időpontokat is mintázzuk (ha a horizonton belül vannak); a félórás rács
  önmagában elmulaszthatná az inhalált koncentráció csúcsát.
  A pillanatképek nem a legközelebbi minta pontatlan kerekítései.
- Dózisidőben N és C folytonos, ezért az összetétel „közvetlenül előtte”
  és „közvetlenül utána” ugyanaz; csak a rejtett tranzitállapot ugrik.
- A fizikai ODE pozitív invariáns: X=0 határon csak nemnegatív mutációs
  beáramlás marad. Ez nem jelenti, hogy a numerikus solver minden belső
  művelete pozitivitásmegőrző.
- Ha egy visszaadott állapot negatívabb, mint −50·atol, vagy nem véges:
  explicit hiba. Kisebb negatív kerekítési maradék nullázása **figyelmeztetéssel
  és exportált számlálóval** történik; nincs elhallgatott negatív populáció.
  atol itt vegyes egységű állapotokra alkalmazott numerikus tolerancia.
  Kis sejtszámok kvantitatív elemzésénél szigorítása/állapotskálázás kell.
- A PD belső solverlépésnél C≥0-t használ, mert tört Hill-kitevő negatív
  koncentrációra nem értelmezhető. Az elfogadott megoldás ettől függetlenül
  átesik a fenti ellenőrzésen.

### Automatikusan ellenőrzött követelmények

| Teszt | Követelmény |
|---|---|
| Logisztikus határeset | Analitikus \(K/[1+(K/N_0-1)e^{-rt}]\), relatív eltérés <2·10⁻⁶ |
| Interakció iránya | A[SA,PA]<0 csak SA-t gátolja a megfelelő kontrollhoz képest |
| Dózisok | 42 aztreonam + 14 azithromycin; 250, illetve 150 mg/beadás; helyes napindexek |
| Gamma PK | Analitikus kéttranzitos koncentrációval egyezés; helyes dózisszuperpozíció |
| PK AUC | β·összdózis/kel, relatív eltérés <10⁻⁵ nagyságrend |
| Mindhárom PK | Független mátrixexponenciális megoldás ellenőrzése |
| Rezisztencia | Gyógyszerspecifikus szelekció, mutáció nullából, fajszintű fluxusmegmaradás |
| Bemenetek | Ismeretlen taxon/szer, hibás dózis, hibás százalék, hiányzó horgony: hiba |
| Others | 65/20/13/2 indulás, változó Others%, mindvégig 100% zárás |
| Kezdeti nulla | Bevándorlás hiányában a nulla faj nem jelenik meg magától |
| Teljes rendszer | 11 taxon, mindhárom szer, 14 nap, véges és nemnegatív kimenet |
| Konvergencia | 0,5→0,25 h maximális lépés és szigorúbb tolerancia: <0,001 log10 CFU/mL eltérés |
| Megőrzés | A korábbi 11 r/K sor és mind a 76 kölcsönhatás pontosan változatlan |
| Export | PNG/CSV/NPZ/JSON létezik, pillanatképek tényleges kért idejűek |

Ezek **szintetikus és implementációs tesztek**, nem betegekhez illeszkedési tesztek.
Nincs jelenleg külső adatokon számított predikciós hiba, profile likelihood,
Morris/Sobol eredmény vagy klinikailag igazolt rezisztenciaidő.

**Az átadott verzión ténylegesen lefutott:** 18/18 unittest sikeres;
mind a négy Python-fájl `py_compile` ellenőrzése sikeres; a notebook 15 cellával,
0 hibakimenettel és 4 beágyazott képpel végrehajtva elmentve.
A helyi dokumentumhivatkozások ellenőrzése nem talált hiányzó fájlt.
A végső Pylance-diagnosztikai lekérdezés az eszköz helyi kapcsolati hibája miatt
nem futott le; teljes végső statikus típusellenőrzési sikert nem állítunk.
A tesztek és a példafutás a jelzett, kismértékű solver-kerekítési korrekciókkal futottak.

## 11. Ábrák és adatexport

A [futtató](szimulacio_futtatas.py) az `eredmenyek` könyvtárba ír:

| Kimenet | Tartalom |
|---|---|
| szimulacios_gorbek.png | Négy panel: log10 CFU/mL, helyi antibiotikum-koncentráció, Δlog10, 100%-os összetétel és dózisjelölések |
| mikrobiom_idopontok.png | Egymás melletti 100%-os oszlopok a kiválasztott napokon; Others szürke |
| rezisztens_frakciok.png | Három gyógyszerhez tartozó fenotípusfrakció taxononként |
| idosor.csv | Nyers CFU, cenzorált log10, százalék, mg/L, gyógyszerspecifikus R-frakció |
| idopontok.csv | A tényleges kért időpontok százalékos összetétele |
| dozisok.csv | Minden beadás ideje és mg/dózis |
| allapotok.npz | Fenotípus-CFU, koncentráció, faj- és bitjelölés, idő, Others; pickle nélkül olvasható |
| futas_beallitasai.json | Tényleges egyenletparaméterek, betegprofil, inicializálás, kezelés, solver és csomagverziók |

A notebook ezen felül `halozat_hipotezis.png` néven az aktuális paraméterekből
készít szűrt hálózati ábrát, és külön `eredmenyek_kontroll` könyvtárba menti a
gyógyszermentes kontrollt.

A koncentrációpanel `symlog` tengelye 0,1 mg/L alatt lineáris; ez mindkét
nagyságrendet és a nullát is láthatóvá teszi. A baktériumok teljes latin
névvel szerepelnek a jelmagyarázatban. A rögzített paletta nem változik
attól, hogy kevesebb taxont választunk. A feliratok a Python-forrásokban
ASCII formájúak; a dokumentáció magyar ékezettel készült.

Az outputkönyvtár azonos nevű eredményeit újrafuttatás felülírja.
Összehasonlításhoz külön `--output` könyvtárat használj. A beállítás-JSON
nyilvántartás, nem automatikusan visszatölthető konfigurációs API.

## 12. Következő validációs lépések, amelyek még nincsenek elvégezve

1. PA és egy szer: gyógyszermentes növekedés, több koncentrációjú time-kill,
   total + szelektív CFU; S/R ráták és PD illesztése.
2. PK külön illesztése ugyanazon mérési mátrixhoz; ne az EC50 és β egymást
   kompenzáló módosításával „magyarázzuk meg” a baktériumgörbét.
3. Kombinációs adatok után explicit szinergia/antagonizmus vizsgálata.
4. Fajok és interakciók fokozatos visszakapcsolása, nem mind a 76 él egyidejű
   becslése egy rövid relatív idősorból.
5. Szintetikus visszanyerhetőség, Morris/Sobol, profile likelihood vagy bootstrap.
6. Betegenként/idősoronként elkülönített külső teszt; abszolút és relatív
   adatok külön megfigyelési hibamodellel, detektálási cenzorálással.

## 13. Források és az ellenőrzés terjedelme

1. Spalding et al. (2018). *Mathematical modelling of the antibiotic-induced
   morphological transition of Pseudomonas aeruginosa.*
   [PLOS Computational Biology, 10.1371/journal.pcbi.1006012](https://doi.org/10.1371/journal.pcbi.1006012).
   Mechanisztikus baktericid tag és tápanyag-/fenotípusmodell; nem az itteni számok forrása.
2. Nashebi, Sari & Kotil (2024). *Mathematical modelling of antibiotic interaction
   on evolution of antibiotic resistance: an analytical approach.*
   [PeerJ 12:e16917](https://pmc.ncbi.nlm.nih.gov/articles/PMC10903357/),
   [DOI](https://doi.org/10.7717/peerj.16917).
   S/R és kombinációs kölcsönhatások koncepcionális alapja.
3. Landersdorfer et al. (2018). *Optimization of a Meropenem-Tobramycin Combination
   Dosage Regimen against Hypermutable and Nonhypermutable Pseudomonas aeruginosa
   via Mechanism-Based Modeling and the Hollow-Fiber Infection Model.*
   [AAC 62:e02055-17](https://pmc.ncbi.nlm.nih.gov/articles/PMC5913933/),
   [DOI](https://doi.org/10.1128/AAC.02055-17).
4. TOBI 300 mg/5 mL SmPC, 5.1. Mechanizmus, CF-köpet bioaktivitása,
   aktivitási és insusceptible lista, adaptív rezisztencia:
   [eMC](https://www.medicines.org.uk/emc/product/262/smpc).
5. Geller et al. (2002). *Pharmacokinetics and bioavailability of aerosolized
   tobramycin in cystic fibrosis.* Chest 122:219–226.
   [10.1378/chest.122.1.219](https://doi.org/10.1378/chest.122.1.219).
   Az absztrakt numerikus PK-adatait ellenőriztük.
6. Cayston EPAR product information, 4.2, 5.1–5.2:
   [EMA PDF](https://www.ema.europa.eu/en/documents/product-information/cayston-epar-product-information_en.pdf);
   [eMC](https://www.medicines.org.uk/emc/product/4456/smpc).
   A PDF 5.2 részében a helyi és szérumadatokat külön ellenőriztük.
7. Wilms et al. (2006). *Pharmacokinetics of azithromycin in plasma, blood,
   polymorphonuclear neutrophils and sputum during long-term therapy in patients
   with cystic fibrosis.* Therapeutic Drug Monitoring.
   [PubMed 16628134](https://pubmed.ncbi.nlm.nih.gov/16628134/),
   [10.1097/01.ftd.0000195617.69721.a5](https://doi.org/10.1097/01.ftd.0000195617.69721.a5).
   A numerikus adatok forrása a publikáció absztraktja, nem saját betegadat.
8. Azithromycin SmPC, különösen 5.2:
   [eMC](https://www.medicines.org.uk/emc/product/6541/smpc).
   Kiegészítő CF-specifikus szisztémás/sejten belüli PK:
   [Absolute bioavailability and intracellular pharmacokinetics… (2005)](https://pmc.ncbi.nlm.nih.gov/articles/PMC1315964/),
   [10.1128/AAC.49.12.5013-5017.2005](https://doi.org/10.1128/AAC.49.12.5013-5017.2005).
9. Gloor et al. (2017). *Microbiome Datasets Are Compositional: And This Is Not
   Optional.* [10.3389/fmicb.2017.02224](https://doi.org/10.3389/fmicb.2017.02224).
10. Vandeputte et al. (2017). *Quantitative microbiome profiling links gut community
    variation to microbial load.* [10.1038/nature24460](https://doi.org/10.1038/nature24460).
    Általános kvantitatív mikrobiom-módszertan, nem CF-specifikus K-forrás.
    A csatolt CFU-dokumentumban is javasolt háttérforrás; a kiadó teljes szövege
    ebben a munkamenetben nem volt hozzáférhető, számszerű paramétert nem vettünk át.
11. *Changes in the lung bacteriome in relation to antipseudomonal therapy in
    children with cystic fibrosis* (2018).
    [10.1007/s12223-017-0562-3](https://doi.org/10.1007/s12223-017-0562-3).
    A cím és DOI bibliográfiailag ellenőrzött; ebből most nem vettünk át számszerű paramétert.
12. Regoes et al. (2004). *Pharmacodynamic Functions: a Multiparameter Approach to
    the Design of Antibiotic Treatment Regimens.*
    [10.1128/AAC.48.10.3670-3676.2004](https://doi.org/10.1128/AAC.48.10.3670-3676.2004).
    A koncentráció–hatás leírás módszertani háttere a fejlesztési tervből.
    A kiadó teljes szövege ebben a munkamenetben nem volt hozzáférhető;
    az implementált Hill/Emax mechanizmust [1–4] is motiválja.

**Forráskritika:** a gyógyszeralkalmazási előírások nem önmagukban gLV/PD
paraméteradatbázisok. Az in vitro kísérlet, a plazma-PK és a CF-köpetből
származó mérés eltérő megfigyelési rendszer; az egyikből a másikba történő
átvitelt minden esetben külön feltevésként kell jelölni.
