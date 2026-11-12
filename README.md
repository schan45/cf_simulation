# CF-mikrobiom szimulátor

Magyar nyelvű, kutatási célú Streamlit felület a TDK mikrobiommodellhez. A modell demonstrációs célú, betegadatokon nem validált, és nem alkalmas klinikai döntéstámogatásra.

## Nyilvános futtatás

Az alkalmazást a GitHubon tárolt forrásból a [Streamlit Community Cloud](https://share.streamlit.io/) tudja elindítani. A GitHub Pages statikus oldalakat szolgál ki, Python/Streamlit alkalmazást nem futtat.

1. Töltsd fel ezt a projektet egy GitHub repositoryba.
2. Nyisd meg a Streamlit Community Cloudot, és jelentkezz be a GitHub-fiókoddal.
3. Válaszd a **Create app** lehetőséget, majd a repositoryt és az `app.py` főfájlt.
4. Telepítés után a Streamlit egy nyilvános `*.streamlit.app` linket ad.

## Helyi futtatás

```powershell
python -m pip install -r requirements.txt
streamlit run app.py
```
