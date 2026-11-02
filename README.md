# Goppa符号復号の実験コード

このリポジトリには、Goppa符号を中心とした符号理論・復号アルゴリズムのC言語による実験実装が含まれています。Berlekamp–Massey法、Euclid法、Patterson法など、異なる復号方法や有限体上の演算を試したコードが複数のディレクトリにあります。

単一のアプリケーションではなく、独立した試作・検証コードの集まりです。ファイルごとに実装状況や前提条件が異なり、複数のソースに `main` 関数が含まれます。

## ディレクトリ

| ディレクトリ | 概要 |
| --- | --- |
| `bm_gf2m/` | GF(2^m)上のBerlekamp–Massey法を用いたGoppa符号復号の実験。 |
| `daikoku/` | Euclid法およびPatterson法によるGoppa符号復号の実験。 |
| `gf2m-peter/` | Goppa符号のPeterson復号器の試作。 |
| `marimo/` | Goppa符号復号と関連する有限体・多項式・行列演算の実験。`vmo.c`ではBerlekamp–Massey系と一般化Patterson復号を試しています。 |
| `marimo/bunka2/`, `marimo/gaal/` | `marimo`内の別バージョンや派生実装。 |

ルート直下にも、これらに関連する試作コードやヘッダがあります。ディレクトリ内のREADMEとMakefileを参照してください。符号パラメータは主に `global.h` または `global-p.h` で定義されています。

## CMakeでのビルド

リポジトリのルートでCMakeを実行すると、登録済みの実験プログラムをまとめてビルドできます。

```sh
cmake -S . -B build
cmake --build build --parallel
```

現在登録しているターゲットは、各コード群の代表的な実行ファイルです。各ディレクトリにあるすべての `.c` ファイルを一括ビルドするものではありません。

| CMakeターゲット | ソース |
| --- | --- |
| `tokyo` | `ymo.c` |
| `bm_gf2m_lyon` | `bm_gf2m/lyon.c` |
| `daikoku_xeuclide` | `daikoku/Xeuclide.c` |
| `gf2m_peter_ptar` | `gf2m-peter/p-tar.c` |
| `marimo_vmo` | `marimo/vmo.c` |

ビルド後は、実行したいターゲットのバイナリを個別に起動します。

```sh
./build/marimo_vmo
```

`daikoku_xeuclide` と `marimo_vmo` はAVX2向けにビルドされます。実行環境のCPUがAVX2に対応している必要があります。実験コードのため、実行に時間がかかったり、大量の診断出力が出たりすることがあります。

従来のMakefileも各ディレクトリに残っています。ルートのMakefileは存在しない `chash-p.h` を依存ファイルに指定しているため、CMakeを使う場合は上記の手順を利用してください。

## 注意

このコードは研究・学習・アルゴリズム検証向けの試作です。実装の完成度や暗号学的安全性は保証されていません。実用の暗号システムや機密データの保護には使用しないでください。
