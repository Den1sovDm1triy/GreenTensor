# SPDX-License-Identifier: MIT
# Scientific scope: scientific research and engineering modeling in classical electrodynamics, antenna theory, microwave devices, and electromagnetic scattering.

"""
ВНИМАНИЕ: этот файл не относится к библиотеке GreenTensor (электродинамика).
Это утилита для анализа HTML соцсетей через Yandex GPT, попавшая сюда по ошибке.

БЕЗОПАСНОСТЬ: ранее здесь был ЗАХАРДКОЖЕН рабочий API-ключ Yandex Cloud
(Api-Key ...) и идентификатор каталога. Ключ удалён из исходника, но он
уже был в репозитории — НЕОБХОДИМО ОТОЗВАТЬ его в Yandex Cloud Console
и выпустить новый. Теперь ключ и folder-id читаются из переменных окружения.
"""
import os
import json
import requests


def analyze_with_gpt(file_path: str):
    api_key = os.environ["YANDEX_API_KEY"]        # экспортируйте перед запуском
    folder_id = os.environ["YANDEX_FOLDER_ID"]

    with open(file_path, "r", encoding="utf-8") as file:
        file_content = file.read()

    url = "https://llm.api.cloud.yandex.net/foundationModels/v1/completion"
    headers = {
        "Content-Type": "application/json",
        "Authorization": f"Api-Key {api_key}",
    }
    data = {
        "modelUri": f"gpt://{folder_id}/yandexgpt/rc",
        "completionOptions": {
            "stream": False,
            "temperature": 0.6,
            "maxTokens": "2000",
        },
        "messages": [
            {
                "role": "system",
                "text": "Ты ассистент, анализируешь контент соцсетей и возвращаешь JSON.",
            },
            {
                "role": "user",
                "text": (
                    "Проанализируй HTML поста и верни JSON с метриками: "
                    f"views, likes, comments, shares.\n\n{file_content}"
                ),
            },
        ],
    }

    response = requests.post(url, headers=headers, json=data)
    if response.status_code == 200:
        result = response.json()
        print("Успешный ответ:")
        print(json.dumps(result, indent=2, ensure_ascii=False))
        return result
    print(f"Ошибка: {response.status_code}")
    print(response.text)
    return None


if __name__ == "__main__":
    analyze_with_gpt("page-raw.txt")
