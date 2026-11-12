# Magyar Streamlit felület

Ez a mappa önálló alkalmazás: az alkalmazás nem importál a szülőkönyvtárból.
A jelenlegi modell, paramétermodul és futtató **másolatai** itt találhatók;
a szülőmappában később módosított értékek nem frissülnek automatikusan itt.

## Indítás

A mappában, megfelelő Python-környezetből:

```powershell
python -m pip install -r requirements.txt
python -m streamlit run app.py --server.address localhost --browser.gatherUsageStats false
```

Ezen a gépen a működő interpreterrel, a projekt gyökeréből:

```powershell
Set-Location .\streamlit_ui
& 'C:\Users\z005awpa\AppData\Local\Programs\Python\Python314\python.exe' -m streamlit run app.py --server.address localhost --browser.gatherUsageStats false
```

Az [inditas.ps1](inditas.ps1) is használható. Ha a `python` rendszeralias nem működik,
állítsd be előtte a `TDK_PYTHON` környezeti változót a Python exe teljes útjára.
A felület alapesetben: http://localhost:8501 . A szerver leállítása: Ctrl+C.
Nem szükséges API-kulcs, felhőszolgáltatás vagy fizetős AI.

## Használat

1. **Páciens:** nem, életkor, BMI metaadat; tüdőkategória és társbetegségek alapján
   demonstrációs immunkonstans. Diabetes és csökkent glükóztolerancia egyszerre nem választható.
2. **Mikrobiom:** legördülő többválasztós lista a 11 taxonhoz, szerkeszthető táblázat.
   - Relatív módban minden fajhoz százalék, legalább egyhez pozitív abszolút
     CFU/mL vagy log10 CFU/mL. A többi abszolút cella maradjon üres.
   - Több horgony esetén 1%-on belüli egyezést ellenőrzünk. Nem egyező adatból
     nem választunk önkényesen „helyeset”.
   - Abszolút módban minden taxon mennyisége kötelező, az Others külön megadható.
   - Relatív módban az Others a százalékok maradéka; a futásban az ebből számított
     abszolút háttér állandó, aránya változó. Az ábrán szürke.
3. **Antibiotikum:** a sor eleji pipa engedélyezi a szert, a napi pipák a kezelési
   napokat jelölik. Napi **összes mg**, beadásszám; a beadások a nap kezdetétől
   egyenletes időközökkel követik egymást. A táblázatos UI egy gyógyszerhez azonos
   napi dózist használ minden kijelölt napon; egyedi napi dózisok a helyi Python
   futtató `DailyDose` API-jával adhatók meg.
4. **Időpontok:** pontosvesszővel elválasztva, például `0; 2,5; 7; 14`.
5. **Futtatás:** a számolás csak a gombra indul. Változtatás után a régi eredmény
   rejtve marad, amíg új futás nem készül.
6. **Eredmény:** görbék, összetételi oszlopok, rezisztens frakciók, kezdeti és végső
   abszolút CFU/mL, változás és opcionális gyógyszermentes kontroll.
   ZIP-ben a tényleges beállítások és minden eredmény letölthető.

**Időhorizont, kiválasztott taxonok, beviteli mód vagy skála váltásakor az érintett
tábla alapértékei újratöltődnek.** Ellenőrizd a cellákat és pipákat futtatás előtt.
Rövid horizontnál a később kezdődő gyógyszer kijelölését is vedd ki, vagy jelölj
neki érvényes napot. Üres kezeléshez mindhárom gyógyszersor kijelölését vedd ki.

## Későbbi paraméterhangolás

A mappán **belüli** [modell_parameterek.py](modell_parameterek.py) szerkesztendő:

- `R_TABLE`: növekedési ráták;
- `K_TABLE`: eltartóképességek;
- `INTERACTIONS`: `(ható taxon, érintett taxon, érték, alsó, felső)`;
- `PK_PARAMETERS`, `PD_PARAMETERS`: gyógyszerparaméterek;
- `IMMUNE_BASE_H`, `LUNG_WEIGHTS`, `DISEASE_WEIGHTS`: immunforgatókönyvek.

Az értékeknek a táblázatok saját alsó–felső határai között kell lenniük.
Ha indokoltan más intervallumra térsz át, azt is módosítani kell.
Mentés után frissítsd a felületet és futtasd újra. Az alkalmazás minden
újrafutáskor újratölti a helyi paramétermodult.
A modell egyenleteit ehhez nem kell módosítani.

## Korlátok és adatkezelés

- **Nincs klinikai siker/nem siker minősítés, „jó/rossz baktérium” automatikus címke
  vagy gyógyszerajánlás.** A végső terhelés csökkenése csak a modell kimenete.
- Az orális/kommensális csoport kijelzése nem jelenti e taxonok garantált előnyösségét.
- Az abszolút eredmény CFU/mL koncentráció, nem a teljes tüdőben lévő CFU.
- A nullából számított százalékos/logaritmikus változás nem értelmezett, üres mezőként jelenik meg.
- Az alkalmazás helyben fut, használati statisztika küldése kikapcsolva.
  Betegazonosítót ne adj meg. A letöltött ZIP-ben a megadott profiladatok szerepelnek.
- Eredmények munkamenetenként külön memóriában; exporthoz saját ideiglenes mappa,
  amely a csomag elkészültekor automatikusan törlődik. Nincs közös eredményfájl-felülírás.
- A helyi szerver nem többfelhasználós éles betegadatkezelő rendszer; ne tedd
  nyilvános hálózatra hitelesítés és megfelelő adatvédelmi tervezés nélkül.
- A másolt [módszertani dokumentum](modell_modszertan.md) a modell alapdokumentuma.
  Egyes történeti notebook/PDF hivatkozásai az eredeti projektben találhatóak,
  nem szükségesek az UI futásához.

## Ellenőrzés

```powershell
python -m unittest test_ui -v
```

A tesztek az adatkonverziót, hibás bemeneteket, napi dózisokat, végső CFU-t,
valamint a Streamlit tényleges felépülését, futtatását és ZIP-exportját ellenőrzik.
A modell eredeti matematikai verifikációja az eredeti projektben külön szerepel.
