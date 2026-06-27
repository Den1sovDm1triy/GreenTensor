# Landing page для greentenor

Одностраничный сайт о библиотеке **GreenTensor** (расчёт рассеяния ЭМ-волн на многослойных сферах методом тензорных функций Грина).

## Структура

```
greentenor_landing/
├── index.html   ← всё в одном файле: HTML + CSS + JS + Three.js-сцена
└── README.md    ← этот файл
```

Страница полностью статична. Все зависимости (Three.js r128, Google Fonts) подгружаются по CDN — никакой сборки и `npm install` не требуется.

## Локальный предпросмотр

Открыть `index.html` двойным кликом — работает, но из-за CORS некоторые браузеры могут не загрузить шрифты при `file://`. Лучше поднять локальный сервер:

```bash
# Python 3 (есть везде)
cd greentenor_landing
python -m http.server 8000
# открыть http://localhost:8000

# Или Node.js
npx serve .
```

## Деплой — варианты по простоте

### Самый простой — Netlify Drop

1. Открыть https://app.netlify.com/drop
2. Перетащить папку `greentenor_landing` в окно браузера
3. Получить URL вида `https://random-name.netlify.app` — сразу работает
4. В настройках сайта → **Domain management** → добавить кастомный домен `greentenor.*` и прописать DNS-записи (Netlify подскажет какие)

### Vercel

```bash
npm i -g vercel
cd greentenor_landing
vercel --prod
```

Кастомный домен подключается через дашборд (Domains → Add).

### Cloudflare Pages

1. Залить папку в любой Git-репозиторий (GitHub / GitLab)
2. Cloudflare Pages → Create project → Connect to Git
3. Build command — оставить пустым, output directory — `/`
4. Кастомный домен — через **Custom domains** в дашборде

### Surge (минимум кликов)

```bash
npm i -g surge
cd greentenor_landing
surge . greentenor.surge.sh
```

### Любой обычный shared-хостинг через FTP/SFTP

Папка `greentenor_landing/` — это и есть готовый «document root». Загрузить её содержимое в `public_html/` (или аналог) на хостинге — и сайт работает.

### GitHub Pages

```bash
cd greentenor_landing
git init
git add .
git commit -m "init: greentenor landing"
git branch -M main
git remote add origin https://github.com/<user>/greentenor.git
git push -u origin main
```

В настройках репозитория → **Pages** → выбрать ветку `main`, папка `/ (root)`.
Подключить домен — через `Custom domain` + добавить файл `CNAME` с содержимым `greentenor.com` (или какой у вас домен).

## Что нужно изменить под себя

| Что | Где |
|---|---|
| Логотип / favicon | в шапке `index.html` (тег `<link rel="icon">` сейчас отсутствует — добавить при необходимости) |
| Email авторов | секция `<footer>`, ссылки `mailto:` |
| Ссылки на GitHub | поиск по тексту `Den1sovDm1triy/GreenTensor` |
| Цветовая палитра | блок `:root` в `<style>`, переменные `--accent`, `--bg`, `--surface` |
| Метаданные SEO | `<meta name="description">`, `<title>` в `<head>` |
| Open Graph / Twitter card | сейчас не добавлены — стоит дописать перед публичным деплоем |

## Размер и производительность

- Один HTML-файл ~28 KB (без сжатия)
- Three.js с CDN ~600 KB (кэшируется браузером)
- Шрифты Google Fonts ~50 KB
- Никаких трекеров, аналитики, кук — добавлять по необходимости

## Лицензия контента

Текст и оформление можно использовать как угодно в рамках вашего проекта. Three.js — MIT.
