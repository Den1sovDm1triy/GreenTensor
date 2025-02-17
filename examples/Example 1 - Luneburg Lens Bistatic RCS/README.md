# GreenTensor Bistatic RCS

## Описание

Этот репозиторий содержит код и результаты моделирования обратного рассеяния (Bistatic RCS) с использованием библиотеки **GreenTensor**. Представленные данные сравниваются с результатами из научной статьи:

**Greenwood и Jin (1999) — "A Novel Efficient Algorithm for Scattering from a Complex BOR Using Mixed Finite Elements and Cylindrical PML"**.

Запустить файл на исполнение можно в [google colab]([url](https://colab.research.google.com/drive/1VpbKq_aoC2bgZgofIQjYkxPGVQFmxZxb#scrollTo=N1dG3kkVsizi)).

## Полученные результаты

На рисунке ниже приведены данные, рассчитанные с использованием GreenTensor:

**Fig. 1 - GreenTensor Bistatic RCS (k₀a=5)**
![GreenTensor Bistatic RCS](fig%201%20-%20green%20tensor%20bistatic%20RCS%20k0a%3D5.png)

Сравнение с литературными данными показывает хорошее соответствие с расчетами из статьи Greenwood и Jin (1999).

## Воспроизведение результатов

Чтобы воспроизвести данные, представленные на **рисунке 11** из статьи Greenwood и Jin (1999), используйте следующий скрипт:

**Luneburg Lens Bistatic RCS**
```bash
python Luneburg\ Lens\ Bistatic\ RCS.py
