# alt_dft

分子・反応データを扱うための Python プロジェクトです。RDKit、RXNMapper、ASE などを利用して、反応マッピング、構造最適化、遷移状態計算、触媒の追加を行います。

## 必要なもの

- [Pixi](https://pixi.sh/)
- `win-64` 環境（現在の `pixi.toml` で対象プラットフォームが固定されています）
- `modal` コマンド（`opt` と `ts` タスクを実行する場合）

Python、RDKit、RXNMapper、ASE、NumPy などの依存関係は Pixi が管理します。

## セットアップ

### pixiの導入
本プロジェクトでは依存関係の管理にpixiを利用しています。
（condaなどの親戚だと思ってください）
pixiをまだ導入していなければ次のコマンドで導入しましょう。
```powershell
winget install --id prefix-dev.pixi
```

### modalの設定
[modal](https://modal.com/)はサーバーレスのGPUクラウドです。
設定をすることで、GPUを借りることができます。
ただし反応系、生成系のマッピングのみを行う場合は必要ありません。

### 初期設定
リポジトリのルート（README.mdがある階層です）で次を実行します。

```powershell
pixi install
```

さらに入力ファイルを置くための`data`フォルダを追加します。
この中で構造ファイルを作成したり最適化したりします。

```powershell
mkdir data
```

### 流れ

まず反応系と生成系のファイルをchemdrawなどで作成します。

ChemDrawならSave AsからMDL Molfile V2000(.mol)として`data`フォルダに保存します。

例えば反応系と生成系のmolファイルを`reactant.mol`, `product.mol`とすれば、`data`フォルダ内で次のコマンドで原子同士を対応させることができます。

```powershell
pixi run map reactant.mol product.mol
```

これにより原子の対応した`.pdb`ファイルが`data`フォルダ内に生成されます。
名前は`reactant_mapped.pdb`, `product_mapped.pdb`です。

これらを構造最適化、TS計算しましょう。
ほどほどにわかってきたら触媒も載せて計算しましょう。

## Pixi コマンド

### プロジェクト固有のタスク

| タスク | 実行内容 |
| --- | --- |
| `map` | `rxn_mapper.py` を実行し、反応系と生成系の原子マッピングを行う |
| `opt` | `modal run opt.py` を実行し、構造最適化を行う |
| `ts` | `modal run ts.py` を実行し、遷移状態関連の計算を行う |
| `putcat` | `pos_cat.py` を実行し、触媒の追加を行う |

実行例：

```powershell
pixi run map <reactant> <product>
pixi run putcat <reactant> <product> <catalyst>
pixi run opt <structure> --charge <charge> --spin <spin>
pixi run ts <reactant> <product> --charge <charge> --spin <spin>
```

### `map` の引数

```text
pixi run map <reactant> <product> [--reactant-out <path>] [--product-out <path>]
```

- `<reactant>`：反応物の入力ファイル
- `<product>`：生成物の入力ファイル
- `--reactant-out`：原子マッピング後の反応物 PDB 出力先（省略時は `data/<入力名>_mapped.pdb`）
- `--product-out`：原子マッピング後の生成物 PDB 出力先（省略時は `data/<入力名>_mapped.pdb`）

入力ファイルは、指定したパスに存在しなければ `data/` 以下から検索されます。

```powershell
pixi run map reactant.mol product.mol
pixi run map data/reactant.mol data/product.mol --reactant-out data/r.pdb --product-out data/p.pdb
```

### `putcat` の引数

```text
pixi run putcat <reactant> <product> <catalyst> [--min-distance <Å>] [--outdir <directory>]
```

- `<reactant>`、`<product>`、`<catalyst>`：PDB または XYZ ファイル。指定パスになければ `data/` 以下から検索されます
- `--min-distance`：基質と触媒の最小原子間距離（Å）。既定値は `8.0`
- `--outdir`：出力先ディレクトリ。既定値は `data`

出力は `<reactant名>_cat.pdb` と `<product名>_cat.pdb` です。反応物と生成物は、原子数と原子順が対応している必要があります。

```powershell
pixi run putcat reactant.pdb product.pdb catalyst.pdb
pixi run putcat reactant.xyz product.xyz catalyst.xyz --min-distance 10 --outdir data/with_catalyst
```

### `opt` の引数

```text
pixi run opt <file> [--charge <整数>] [--spin <整数>] [<file> [--charge <整数>] [--spin <整数>] ...]
```

入力は PDB または XYZ で、指定パスになければ `data/` 以下から検索されます。1 回の実行で複数ファイルを指定できます。各ファイルの直後にオプションを置く形式です。

- `--charge`：分子の電荷。既定値は `0`
- `--spin`：スピン多重度。既定値は `1`

最適化結果は `data/optimized/<入力名>_opt.pdb` に保存されます。

```powershell
pixi run opt molecule.xyz
pixi run opt molecule.xyz --charge -1 --spin 2
pixi run opt a.pdb --charge 0 b.xyz --charge 1 --spin 2
```

### `ts` の引数

```text
pixi run ts <reactant> <product> [--charge <整数>] [--spin <整数>] [<reactant> <product> ...]
```

反応物と生成物を 1 組として指定します。PDB または XYZ を使用でき、指定パスになければ `data/` 以下から検索されます。複数反応を続けて指定できます。

- `--charge`：その反応の電荷。既定値は `0`
- `--spin`：その反応のスピン多重度。既定値は `1`

遷移状態候補は `data/TS_candidates/<反応物名>_TS.log` に保存されます。オプションは対応する反応物・生成物の直後に置いてください。

```powershell
pixi run ts reactant.pdb product.pdb
pixi run ts reactant.xyz product.xyz --charge 0 --spin 1
pixi run ts r1.pdb p1.pdb --charge -1 r2.pdb p2.pdb --spin 2
```

`opt` と `ts` は Pixi の環境から Modal CLI を呼び出します。Modal の認証や実行環境が必要な場合は、Modal 側の設定を先に済ませてください。

## 主なファイル

- `rxn_mapper.py`：反応の原子マッピング
- `pos_cat.py`：触媒候補の処理
- `molecule_io.py`：分子構造の入出力
- `atom_mapping.py`：原子マッピング関連の処理
- `gaussian.py`：Gaussian 入力などの処理
- `opt.py`：構造最適化の Modal タスク
- `ts.py`：遷移状態計算の Modal タスク
- `data/`：入力データ・生成データの保存場所
- `pixi.toml`：依存関係とタスクの定義
- `pixi.lock`：依存関係の固定情報
