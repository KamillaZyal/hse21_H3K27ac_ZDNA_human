# hse21_H3K27ac_ZDNA_human
# Анализ пиков гистоновой метки
## Первичная обработка файлов
- Загрузим 2 .bed файла ChIP-seq экспериментов из ENCODE (пики гистоновой метки):
 ```
 wget https://www.encodeproject.org/files/ENCFF389RXK/@@download/ENCFF389RXK.bed.gz
 wget https://www.encodeproject.org/files/ENCFF926NKP/@@download/ENCFF926NKP.bed.gz
 ```
- Распакуем скачанные архивы, оставив необходимые первые 5 столбцов .bed файлов:
```
zcat ENCFF832EOL.bed.gz  |  cut -f1-5 > H3K4me3_A549.ENCFF832EOL.hg38.bed (аналогично для второго файла)
```
- Переведем координаты ChIP-seq пиков к версии генома hg19 с помощью liftOver:
     - Загрузим необходимый файл для конверации hg38 в hg19:
    ```
    wget https://hgdownload.cse.ucsc.edu/goldenpath/hg38/liftOver/hg38ToHg19.over.chain.gz
    ```
    - Запустим liftOver на файлах .bed:
    >Используем команду в формате: *liftOver input.bed hg18ToHg19.over.chain.gz output.bed unlifted.bed*,
    >>где input.bed - файл с исходными данными, hg18ToHg19.over.chain.gz - вспомогательный файл для конвертации,
    >>output.bed - файл, содержащий конвертированные координаты, unlifted.bed будет содержать все позиции генома, которые нельзя перевести в нужную вурсию генома
     ```
    liftOver   H3K4me3_A549.ENCFF832EOL.hg38.bed   hg38ToHg19.over.chain.gz   H3K4me3_A549.ENCFF832EOL.hg19.bed   H3K4me3_A549.ENCFF832EOL.unmapped.bed (аналогично для второго     файла)
    ```
- Построим гистограмму длин участков для каждого эксперимента до и после конвертации к нужной версии генома:
  > Для обработки данных используем скрипт */src/filter_peaks.R*
    - ***Гистограммы для версии генома hg38***
    ![Alt-текст](https://github.com/KamillaZyal/hse21_H3K27ac_ZDNA_human/blob/main/images/filter_peaks.H3K27ac_A549.ENCFF389RXK.hg38.init.hist.png) 
    ![Alt-текст](https://github.com/KamillaZyal/hse21_H3K27ac_ZDNA_human/blob/main/images/filter_peaks.H3K27ac_A549.ENCFF926NKP.hg38.init.hist.png)
    - ***Гистограммы для версии генома hg19***
    ![Alt-текст](https://github.com/KamillaZyal/hse21_H3K27ac_ZDNA_human/blob/main/images/filter_peaks.H3K27ac_A549.ENCFF389RXK.hg19.init.hist.png) 
    ![Alt-текст](https://github.com/KamillaZyal/hse21_H3K27ac_ZDNA_human/blob/main/images/filter_peaks.H3K27ac_A549.ENCFF926NKP.hg19.init.hist.png)
- Проведем фильтрацию пиков по длине:
    - Рассмотрим и сравним самые длинные пики для разных версий геномов:
    > |  | Версия генома hg38 | Версия генома hg19 |
      |----------------|:---------:|:---------:|
      | ENCFF926NKP| 118990 |  chrom     start       end      name score   len |
                            |1 chr20  45927842  45992587  Peak_191  1000 64745 |
                            |2 chr16  87854991  87906505 Peak_4576  1000 51514 |
                            |3  chr3  53260249  53309750 Peak_1352  1000 49501 |
                            |4 chr15  99392468  99441956 Peak_1710  1000 49488 |
                            |5 chr22  30615880  30658937  Peak_660  1000 43057 |
                            |6  chr5 172321108 172362606   Peak_86  1000 41498 |
      | ENCFF389RXK| 118767 | 115997 | 
  >Установим порог 5kb (не длиннее 5000 пар оснований)
  >>Для обработки данных используем скрипт */src/filter_peaks.R - раздел Remove long peaks*
    - ***Гистограммы для filter_peaks***
    ![Alt-текст](https://github.com/KamillaZyal/hse21_H3K27ac_ZDNA_human/blob/main/images/filter_peaks.H3K27ac_A549.ENCFF389RXK.hg19.filtered.hist.png) 
    ![Alt-текст](https://github.com/KamillaZyal/hse21_H3K27ac_ZDNA_human/blob/main/images/filter_peaks.H3K27ac_A549.ENCFF926NKP.hg19.filtered.hist.png)
- Статистика по количеству пиков:
  |  | ENCFF926NKP | ENCFF389RXK |
  |----------------|:---------:|:---------:|
  | Для версии генома ***hg38*** | 118990 | 116220 |
  | Для версии генома ***hg19*** | 118767 | 115997 | 
  | Для ***filter_peaks*** | 115115 | 112077 |  

    
