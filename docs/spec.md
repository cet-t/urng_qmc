# urng\_qmc\_core / urng\_qmc 仕様書

Oct 14, 2026 · @cet（同日改訂: レビュー反映。変更点と根拠は第 11 章）

## 1. 概要

準乱数列（低食い違い列）を 2 つの `no_std` クレートに分けて提供する。`urng_qmc_core` はトレイト（`QuasiRng<D>`、`Seek`、`DigitalNet<D>`、`QmcWord`）と、トレイトだけで書ける薄い接着剤（`AsRng`）を持ち、`urng_qmc` がアルゴリズム（Weyl、R\_d、Sobol など）とスクランブルを実装する。1 次元の列は `AsRng` で包むと `urng_core::Rng` として使える。

分ける理由は 2 つ。独自の列を作る利用者がアルゴリズム群とその方向数表に依存せずに済むこと、そしてトレイトだけを早く安定させ、アルゴリズム側は破壊的変更を気軽に入れられること。

**目的**

- `urng_core` の利用者が、生成器を作る箇所を差し替えるだけで擬似乱数から準乱数に移行できること（1 サンプルにつき 1 値を引くコードに限る。理由は 4.4）
- 多次元列を `nextu()` の連続呼び出しで誤用できない API にすること
- `urng_core` の擬似乱数をスクランブル源として使い、RQMC（ランダム化準モンテカルロ）を組めること

**スコープ（v0.1〜v0.3）**

- 列: Weyl / R\_d（Roberts）、van der Corput（基底2）、Sobol（Joe–Kuo）、Halton
- スクランブル: Cranley–Patterson 回転、デジタルシフト（XOR）
- 任意インデックスへのジャンプ（`Seek`）

**非目標**

- 暗号用途（`urng_core` と同じく非対応）
- `alloc` を要する機能、実行時に次元数が決まる列
- 正規分布など \[0,1) 以外への分布変換（利用者側の責務。開区間 (0,1) への写像は 7.1 の将来候補）
- `urng_core` 本体の変更（本クレートは既存 API のみで成立させる。`Word` に足りない演算は `QmcWord` で補う）

## 2. 前提: urng\_core の変換は低食い違い性を壊さない

準乱数の出力は「\[0,1) の固定小数点数」として `u32` / `u64` に載せる。`urng_core` の変換はすべて**上位ビット基準で単調非減少**なので、列の均一性がそのまま出力側に残る。

| 関数 | 実装 | QMC への影響 |
| --- | --- | --- |
| `Word::to_float` | 上位 23 / 52 bit を仮数部に詰めて `- 1.0` | 下位ビットの切り捨てのみ。解像度は f32 で 2^-23、f64 で 2^-52（u64 は urng\_core 0.1.2 以降で正しい。下記） |
| `Word::randi` | 倍幅 multiply-shift | 単調な写像なので層別構造を保つ。`%` 方式なら壊れていた |
| `Word::randf` | `to_float` の線形変換 + `next_down` | 線形なので保たれる |

**要件**: 各列は「値の上位ビットほど均一性が高い」表現で出力すること。Sobol や van der Corput のように下位ビットが 0 になる列も、`to_float` は上位から使うので問題ない。

**urng\_core 0.1.1 の不具合（0.1.2 で修正済み）**: `u64::to_float` が `x >> 11` で仮数部を作るため、最上位ビットが指数部に落ちて消える（f64 の符号 + 指数は 12 bit なので `x >> 12` が正しい）。結果は x / 2^64 ではなく (x mod 2^63) / 2^63 になり、`0x8000…` が 0.5 ではなく 0.0、`0x4000…` が 0.25 ではなく 0.5 になる。`randf` も同じ変換を使う。擬似乱数では一様性が保たれるので表に出にくいが、準乱数では点が x → 2x mod 1 に折り畳まれ、2^m 点の各座標が 2 重に重なる。u32 は正しい。

- 0.1.2 で `(64 bias) => 12` に修正された。両クレートは `urng_core = "0.1.2"` を要求する
- 再発防止として、`urng_qmc` のテストで `next_point_f` が u64 の上位 52 bit / u32 の上位 23 bit を使うことを確かめる。discrepancy などの品質テストは `to_float` に依存しないよう、浮動小数を自前で計算する

**注意**: `to_float` は f32 で 23 bit しか使わない。`u32` の Sobol を 2^23 点より多く f32 化すると、隣接点が同じ f32 に重なり始める。大量に使う場合は `u64` 版を推奨とドキュメントに明記する。

## 3. クレート構成

トレイトとアルゴリズムを別クレートにし、`urng_core` と合わせて 3 段の依存にする。両クレートとも `#![no_std]`、`alloc` 不要、`#![forbid(unsafe_code)]`。

```
urng_core  ←  urng_qmc_core  ←  urng_qmc
     ↑______________________________|
```

`urng_qmc_core` は依存を `urng_core` だけに保ち、列の実装を一切持たない。独自の列を作る利用者は `urng_qmc` を入れずに済む。

### 3.1 urng\_qmc\_core

```toml
[package]
name = "urng_qmc_core"
version = "0.1.0"
edition = "2024"
rust-version = "1.86"
license = "MIT OR Apache-2.0"
keywords = ["no_std", "qmc", "low-discrepancy", "traits"]

[dependencies]
urng_core = "0.1.2"
```

### 3.2 urng\_qmc

```toml
[package]
name = "urng_qmc"
version = "0.1.0"
edition = "2024"
rust-version = "1.86"
license = "MIT OR Apache-2.0"
keywords = ["no_std", "qmc", "sobol", "low-discrepancy"]

[dependencies]
urng_core = "0.1.2"
urng_qmc_core = { version = "0.1", path = "../urng_qmc_core" }

[features]
default = []
sobol-21201 = []   # SobolWide32/64（Joe–Kuo 全 21201 次元）を足す（5.4 参照）
```

- `path` はローカル開発用。`cargo publish` 時は `path` が外れて `version` だけが残るので、`urng_qmc_core` を先に公開する
- Sobol は v0.1 から feature なしで 1〜1024 次元を扱う。パラメータ表はコンパイル時の `const` 評価にしか使わず**バイナリには入らない**ため、次元上限を feature で絞る利点はコンパイル時間だけで、1024 次元までは無視できる（5.4）

### 3.3 リポジトリ構成

クレートごとに別リポジトリ（`cet-t/urng_qmc_core`、`cet-t/urng_qmc`）とし、ローカルでは 2 つを兄弟ディレクトリに clone する。CI の `urng_qmc` 側は両リポジトリを兄弟ディレクトリに checkout してから実行する（両方 public なのでトークン不要）。

```
urng_qmc_core/
├─ Cargo.toml
└─ src/
   ├─ lib.rs               // 再エクスポート、Float<W> 型エイリアス
   ├─ word.rs              // QmcWord（sealed）
   ├─ quasi.rs             // QuasiRng<D>
   ├─ seek.rs              // Seek
   ├─ net.rs               // DigitalNet<D>
   └─ adapter.rs           // AsRng<Q>
urng_qmc/
├─ Cargo.toml
├─ NOTICE                  // Joe–Kuo 方向数の著作権表示
├─ docs/spec.md            // 本書
├─ src/
│  ├─ lib.rs               // pub use urng_qmc_core::*; と各型の再エクスポート
│  ├─ weyl.rs              // Weyl32/64
│  ├─ rd/
│  │  ├─ mod.rs            // Rd32<D>/Rd64<D>
│  │  └─ table.rs          // R_d 定数（生成物、手で編集しない）
│  ├─ vdc.rs               // VanDerCorput32/64
│  ├─ sobol/
│  │  ├─ mod.rs            // Sobol32<D>, Sobol64<D>
│  │  ├─ wide.rs           // SobolWide32<D>, SobolWide64<D>（sobol-21201）
│  │  ├─ table.rs          // Joe–Kuo パラメータ 1024 次元（生成物、手で編集しない）
│  │  └─ table_21201.rs    // 同 全 21201 次元（生成物、sobol-21201 のときだけ使う）
│  ├─ halton.rs            // Halton32<D>, Halton64<D>（v0.3）
│  └─ scramble.rs          // Shifted, DigitalShifted（v0.2）
├─ tests/                  // 結合テスト（std 可）。既知値・層別・discrepancy・収束
│  └─ data/                // 期待値（生成物）
└─ tools/
   ├─ gen_rd.py            // R_d 定数を 80 桁の decimal で生成
   ├─ gen_sobol.py         // new-joe-kuo-6.21201（または scipy 同梱の同データ）→ table.rs
   └─ gen_vectors.py       // scipy.stats.qmc.Sobol から期待値を生成
```

### 3.4 切り分けの方針

| 置き場 | 入れるもの | 理由 |
| --- | --- | --- |
| `urng_qmc_core` | `QuasiRng<D>`、`Seek`、`DigitalNet<D>`、`QmcWord`、`Float<W>`、`AsRng<Q>` | 外部の列が実装・利用する契約と、契約だけで書ける接着剤。`DigitalNet` も契約なのでここに置き、独自のデジタルネットも `DigitalShifted` を使えるようにする。`AsRng` は外部の 1 次元列も `Rng` 化できるようにここに置く |
| `urng_qmc` | 全列、`Shifted`、`DigitalShifted`、定数表 | スクランブルもトレイトだけで書けるが、列の変換ロジックは持たないという core 側の原則を優先する |

- **再エクスポート**: `urng_qmc` は `pub use urng_qmc_core::*;` する。利用者は通常 `urng_qmc` だけを依存に書けばよい
- **バージョン**: `urng_qmc_core` は早めに 1.0 を目指し、破壊的変更を避ける。core のメジャーが上がると、古い core に依存する第三者の列と新しい `urng_qmc` のトレイトが別物になるため
- **公開依存**: `urng_core::Word` / `Rng` が core の公開 API に現れるので、`urng_core` は core の公開依存になる。`urng_core` 0.2 が出ると core も破壊的変更になるため、core の 1.0 は `urng_core` の 1.0 を前提にする
- **sealed**: `QmcWord` は sealed にし、後からメソッドや上位トレイトを足しても破壊的変更にならないようにする
- **孤児ルール**: `impl Rng for AsRng<Q>` は `AsRng` が core 内の型なので書ける。`impl<Q: QuasiRng<1>> Rng for Q` のような包括実装は、`Rng` も `Q` も core の外のものなので孤児ルール（E0210）で書けない
- MSRV は両クレートとも `urng_core` に合わせて 1.86。定数生成は `const fn` で行い、`build.rs` は使わない
- `pastey` への直接依存は持たない。u32/u64 の展開は内部 `macro_rules!` で済ませる

## 4. 公開 API

4.1 のトレイトは `urng_qmc_core` で定義し、`urng_qmc` が再エクスポートする。4.2 の型はすべて `urng_qmc` に置く。`Word` が sealed なので、各列はジェネリックな `W` ではなく `…32` / `…64` の具象型をマクロで生成する。

### 4.1 トレイト

```rust
// urng_qmc_core
use urng_core::{Rng, Word};

pub type Float<W> = <W as Word>::Float;

/// QMC で使うワード。urng_core::Word に、ジェネリックなラッパーが要る演算を足したもの。
/// sealed。u32 と u64 だけが実装する。
pub trait QmcWord:
    Word + Copy + Eq + Ord + Hash + Debug + Default + Send + Sync + 'static
    + BitXor<Output = Self> + BitXorAssign
{
    const BITS: u32;
    fn wrapping_add(self, rhs: Self) -> Self;
    /// n mod 2^BITS。
    fn from_index(n: u64) -> Self;
    fn to_u64(self) -> u64;
    /// 任意幅の Rng から 1 ワード引く（32→64 は 2 回引いて連結、64→32 は上位を使う）。
    fn from_rng<R: Rng + ?Sized>(rng: &mut R) -> Self where R::Word: QmcWord;
}

/// D 次元の低食い違い列。各成分は [0,1) の固定小数点数。
pub trait QuasiRng<const D: usize> {
    type Word: QmcWord;

    /// 次の点を返し、インデックスを 1 進める。
    #[must_use]
    fn next_point(&mut self) -> [Self::Word; D];

    /// next_point と同じ点を out に書く。D が大きい列はスタック上のコピーを避けるために上書きする。
    #[inline]
    fn fill_point(&mut self, out: &mut [Self::Word; D]) { *out = self.next_point(); }

    /// 各成分を [0,1) の浮動小数にして返す。
    #[must_use]
    #[inline]
    fn next_point_f(&mut self) -> [Float<Self::Word>; D] {
        self.next_point().map(Word::to_float)
    }

    /// 各成分を [min[i], max[i]) に写して返す。
    #[must_use]
    #[inline]
    fn next_point_in(
        &mut self,
        min: [Float<Self::Word>; D],
        max: [Float<Self::Word>; D],
    ) -> [Float<Self::Word>; D] {
        let p = self.next_point();
        core::array::from_fn(|i| p[i].randf(min[i], max[i]))
    }
}

/// 任意インデックスへのジャンプ。並列化や再開に使う。
pub trait Seek {
    /// 次の next_point / nextu が index 番目の点を返す状態にする。index は mod 2^B に丸める。
    fn seek(&mut self, index: u64);
    /// 次に返す点のインデックス。常に [0, 2^B)。
    fn index(&self) -> u64;
}

/// 1 次元の列を urng_core::Rng として使うための明示的なアダプタ。
pub struct AsRng<Q>(/* private */);
impl<Q> AsRng<Q> {
    pub const fn new(inner: Q) -> Self;
    pub fn into_inner(self) -> Q;
    pub const fn get_ref(&self) -> &Q;
    pub fn get_mut(&mut self) -> &mut Q;
}
impl<Q: QuasiRng<1>> Rng for AsRng<Q> { type Word = Q::Word; /* nextu = next_point()[0] */ }
impl<Q: Seek> Seek for AsRng<Q> { /* 委譲 */ }
```

- 整数範囲への写像（`randi` 相当）はトレイトに入れない。利用者が `p[i].randi(..)` を呼べば済む
- ジェネリックに浮動小数を扱いたい利用者は `Float<Q::Word>: Into<f64>` を境界に足せばよい（f32 / f64 とも満たす）
- `fill_point` は既定実装つきなので、後から最適化しても実装側の破壊的変更にならない

### 4.2 型一覧

| 型 | D | `Seek` | `DigitalNet` | 版 |
| --- | --- | --- | --- | --- |
| `Weyl32`, `Weyl64` | 1 | あり | なし | v0.1 |
| `Rd32<D>`, `Rd64<D>` | 1〜16 | あり | なし | v0.1 |
| `VanDerCorput32`, `VanDerCorput64` | 1 | あり | あり | v0.1 |
| `Sobol32<D>`, `Sobol64<D>` | 1〜1024 | あり | あり | v0.1 |
| `SobolWide32<D>`, `SobolWide64<D>` | 1〜21201（feature `sobol-21201`） | あり | あり | v0.3 |
| `Shifted<Q, D>` | Q に従う | Q が実装すれば | なし | v0.2 |
| `DigitalShifted<Q, D>` | Q に従う | Q が実装すれば | あり | v0.2 |
| `Halton32<D>`, `Halton64<D>` | 1〜32 | あり | なし | v0.3 |

どの型も `Rng` を直接は実装しない。`D = 1` の列は `AsRng::new(q)` で `Rng` になる。

### 4.3 使用例

```rust
use urng_qmc::{QuasiRng, Sobol64};

let mut q = Sobol64::<2>::new();
let mut acc = 0.0f64;
for _ in 0..(1 << 16) {
    let [x, y] = q.next_point_f();
    if x * x + y * y < 1.0 { acc += 1.0; }
}
let pi = 4.0 * acc / (1 << 16) as f64;
```

```rust
use urng_core::Rng;
use urng_qmc::{AsRng, Sobol64};

// 1 サンプルにつき 1 値だけ引く既存コードなら、生成器を作る箇所だけ差し替えればよい
let mut rng = AsRng::new(Sobol64::<1>::new());
let u: f64 = rng.nextf();
```

### 4.4 1 次元列を `Rng` として直接実装しない理由

`Rng` を受け取る既存コードの多くは 1 サンプルで複数の値を引く（2 次元の点、Box–Muller、棄却法、シャッフルなど）。1 次元の低食い違い列の連続値は強く相関しているので、これらは黙って壊れる。

実測（N = 2^16 サンプル、連続する 2 値を (x, y) として 4·P(x²+y² < 1) を計算）:

| 生成器 | 推定値 |
| --- | --- |
| van der Corput / Sobol\<1\> の連続値 | 3.2915 |
| Weyl（黄金比）の連続値 | 3.2674 |
| `Sobol<2>` の点 | 3.1418 |

型を差し替えるだけで `Rng` になると、この誤用がコンパイルを通る。`AsRng` で明示的に包ませ、そのドキュメントに上記の制約を書く。あわせて、同じ型が `Rng` と `QuasiRng<1>` の両方を実装すると、ジェネリックな `Q: QuasiRng<1> + Rng` で `Q::Word` が曖昧になる（E0221）問題も避けられる。

## 5. 各列の仕様

全列で共通の約束: 最初の点のインデックスは 0。`next_point` は「現在の点を返してから進める」（後置インクリメント）。B はワード幅（32 / 64）。インデックスは常に \[0, 2^B) に保ち、加算はラップさせる（`overflow-checks = true` の release でも panic しない）。

全型 `Default` を実装する（`new()` / `golden()` と同じ）。

### 5.1 Weyl（`Weyl32` / `Weyl64`）

x\_n = (s₀ + n·α) mod 2^B。

- 状態: `state: W`, `alpha: W`, `seed: W`（= s₀）, `index: u64`
- `new(alpha)`: 常に `alpha | 1` に補正する（debug / release で挙動を変えない）。奇数でないと周期が 2^B 未満になる。どんな奇数でも低食い違い性が出るわけではない（α = 1 なら x\_n = n）ので、通常は `golden()` を使うようドキュメントに書く
- `golden()`: α = 黄金比の逆数（5.2 の表の D = 1）
- `with_seed(self, s₀) -> Self`: 開始点をずらす。既定は s₀ = 0。現在のインデックスは保つ
- `seek(n)`: `state = seed.wrapping_add(alpha.wrapping_mul(n as W))`、O(1)

### 5.2 R\_d（`Rd32<D>` / `Rd64<D>`）

Roberts 列。次元 j（1〜D）の成分はα\_j を使う Weyl 列で、φ\_D は x^(D+1) = x + 1 の正の実根。

```latex
\alpha_j = \left\lfloor \varphi_D^{-j} \cdot 2^B \right\rfloor \;\mathrm{OR}\; 1
```

- 定数は `tools/gen_rd.py` が 80 桁の `decimal` 演算で u64 を生成し、D = 1〜16 を `const` 表として埋め込む。u32 版は u64 の上位 32 bit に OR 1 したもの（⌊α·2^32⌋ = ⌊α·2^64⌋ >> 32 なので一致する）
- 丸めではなく「切り捨ててから最下位を 1」にする。丸めると D = 1 の u64 が `…7C16`、D = 3 の j = 1 が `…ED04`（ともに偶数）になり、周期が半分になるため
- 状態: `state: [W; D]`, `seed: [W; D]`, `index: u64`。`with_seed(self, [W; D])` で開始点をずらす
- D が 1〜16 の範囲外なら `const { assert!(..) }` でコンパイルエラー

| D | j | u64 | u32 |
| --- | --- | --- | --- |
| 1 | 1 | `0x9E3779B97F4A7C15` | `0x9E3779B9` |
| 2 | 1 | `0xC13FA9A902A6328F` | `0xC13FA9A9` |
| 2 | 2 | `0x91E10DA5C79E7B1D` | `0x91E10DA5` |
| 3 | 1 | `0xD1B54A32D192ED03` | `0xD1B54A33` |
| 3 | 2 | `0xABC98388FB8FAC03` | `0xABC98389` |
| 3 | 3 | `0x8CB92BA72F3D8DD7` | `0x8CB92BA7` |

この表は生成スクリプトのテスト期待値にも使う（mpmath 80 桁で全値一致を確認済み）。

### 5.3 van der Corput（`VanDerCorput32` / `VanDerCorput64`）

基底 2 の径向逆関数。x\_n = `(n as W).reverse_bits()`。

- 状態は `index: u64` のみ。`seek` は代入だけ
- 実質的に Sobol の 1 次元目と同じ点集合（順序は Gray code 分異なる）

### 5.4 Sobol（`Sobol32<D>` / `Sobol64<D>`）

Joe–Kuo `new-joe-kuo-6.21201` の方向数を使う。Gray code 更新で 1 点あたり O(D)。

**方向数の構築**（次元ごと、k = 0〜B−1）

- 次元 1: v\_k = 1 << (B−1−k)
- 次元 2 以降: 表から次数 s、係数 a（s−1 bit）、初期値 m₁〜m\_s を読む
  - k < s: v\_k = m\_(k+1) << (B−1−k)
  - k ≥ s: 下の漸化式。a\_i = (a >> (s−1−i)) & 1

```latex
v_k = v_{k-s} \oplus (v_{k-s} \gg s) \oplus \bigoplus_{i=1}^{s-1} a_i \, v_{k-i}
```

**メモリ配置**: 方向数は `impl<const D: usize> Sobol64<D> { const V: &'static [[u64; D]; 64] = &build::<D>(); }` として `const fn` で生成し、rodata に置く。

- `const V: [[u64; 64]; D]`（値の const）にしてはいけない。const は使うたびに一時値として実体化されるため、debug ビルドでは添字アクセスのたびに表全体がスタックへコピーされる（D = 1024 で 512 KB のフレームとループ内 `memcpy` を確認）。`&'static` にすると 120 B のフレームで済む
- 添字は `V[k][d]`（ビット位置が外側）。更新で触る V\[c\] が連続領域になり、ベクトル化しやすい
- インスタンスの状態は `x: [W; D]` と `index: u64` だけになる
- パラメータ表は `const` 評価でしか使わないのでバイナリには入らない。バイナリに入るのは使った D の V だけ（D·B·B/8 バイト。u64・D = 1024 で 512 KB）

**更新**

```rust
fn next_point(&mut self) -> [u64; D] {
    let out = self.x;
    let c = self.index.trailing_ones() as usize; // 最下位の 0 bit の位置
    if c == 64 {                                  // 2^B 点を使い切った
        self.x = [0; D];
        self.index = 0;
    } else {
        for d in 0..D { self.x[d] ^= Self::V[c][d]; }
        self.index += 1;
    }
    out
}
```

**seek(n)**: n を mod 2^B に丸め、g = n ^ (n >> 1) とし、g の立っているビット k について V\[k\] を XOR する。O(D·B)。

**参照値**（D = 2、先頭 8 点、`to_float` 後）: (0, 0), (0.5, 0.5), (0.75, 0.25), (0.25, 0.75), (0.375, 0.375), (0.875, 0.875), (0.625, 0.125), (0.125, 0.625)（`scipy.stats.qmc.Sobol(2, scramble=False)` と一致を確認済み）

**次元の上限**

- `Sobol32/64<D>` は表に 1024 次元（次数 s ≤ 13）を持ち、feature なしで `D ≤ 1024` を受け付ける
- D ごとに表全体を 1 回の `const` 評価で作る方式は、数千次元で deny 既定の `long_running_const_eval` リントに当たる（u64 で D = 2048 は通り、4096 で止まる）。しかもこのリントの許可レベルは、評価を実体化したクレート（利用者側）のルートで決まる（評価中のフレームが他クレートのものだと、レベルを判定する場所がクレートルートに落ちるため）。V や `build` に `#[allow]` を付けても効かず、利用者に `#![allow(long_running_const_eval)]` を書かせることになる
- そこで `sobol-21201` feature では別の型 `SobolWide32/64<D>` を足す。feature は型を足すだけなので加算的で、`Sobol32/64` の挙動は変わらない
  - 方向数を 256 次元ずつのチャンクに分け、チャンクごとに別の関連 const（`Chunk::<C>::V64`）にする。1 回の評価がリントの上限に届かない
  - 型ごとのチャンク表 `CHUNKS` は `SobolWide<D>` の関連 const で、`ceil(D/256)` 個のチャンクだけを評価し、残りはゼロのチャンクで埋める。チャンクは関数呼び出し越しに参照する（初期化式に const を直接書くと、通らない分岐の分まで評価されるため）
  - ジェネリックでない const は定義側のクレートで必ず評価されるので、チャンク表を自由な const にしてはいけない（試作では feature を有効にしただけで `urng_qmc` のビルドが約 200 秒延びた）
  - 実測のビルド時間: feature を有効にして `SobolWide` を使わない場合は増えない。`SobolWide64<3000>` は 8 秒、全 21201 次元を両方の幅で使うテストは約 2 分。バイナリに入るのは使うチャンクだけ
  - 更新はチャンクごとに V\[c\] の行を XOR する。`D ≤ 1024` では `Sobol32/64` と同じ点を返す（テストで確認）
- D = 21201 では `next_point` の戻り値が 170 KB になるので、`fill_point` を勧める

### 5.5 Halton（`Halton32<D>` / `Halton64<D>`、v0.3）

次元 j の基底は j 番目の素数（2, 3, 5, …, 131）。径向逆関数 φ\_b(n) を固定小数点で返す。

- 両幅とも、桁ごとに d\_k · ⌊(2^128 − 1) / b^k⌋ を u128 で加算し、上位 B bit を返す。誤差は 2^-118 程度で、出力は真値の切り捨てとほぼ常に一致する（基底 3, 5, 7, 131 で各 2 万点、誤差 0 を確認）
  - 「u64 で ⌊2^64 / b^k⌋ を加算」する方式は誤差が (b−1)·K ulp（K は桁数）まで積もる。基底 131 では上限 1300 ulp、実測でも 559 ulp あり、`to_float` の解像度（4096 ulp）に対して余裕がない
  - u128 の除算を毎点行う方式（`(num << 32) / den`）は、thumbv7em などでソフトウェア除算になり遅い
- インデックスは \[0, 2^B) でラップする（Halton 自体に周期はないが、他の列と `Seek` の約束を揃える）
- 基底が大きい次元は先頭の点が強く相関する。D > 8 では Sobol を勧めるとドキュメントに書く

## 6. スクランブル（v0.2）

スクランブルは「任意の `urng_core::Rng` から乱数を D 個引いて、列の全点に同じ変換をかける」ラッパーとして実装する。異なるシードで M 個の独立な推定を作れば、そのばらつきから誤差を推定できる（RQMC）。

```rust
pub struct Shifted<Q, const D: usize> where Q: QuasiRng<D> { inner: Q, shift: [Q::Word; D] }
pub struct DigitalShifted<Q, const D: usize> where Q: DigitalNet<D> { inner: Q, mask: [Q::Word; D] }

impl<Q: QuasiRng<D>, const D: usize> Shifted<Q, D> {
    pub fn new<R: Rng + ?Sized>(inner: Q, rng: &mut R) -> Self where R::Word: QmcWord;
    pub const fn with_shift(inner: Q, shift: [Q::Word; D]) -> Self;
    pub const fn shift(&self) -> &[Q::Word; D];
    pub const fn get_ref(&self) -> &Q;
    pub fn into_inner(self) -> Q;
}
// DigitalShifted も同じ形（new / with_mask / mask / get_ref / into_inner）

/// デジタルネット（基底 2）である列のマーカー。定義は urng_qmc_core。
/// urng_qmc では Sobol と van der Corput が実装し、外部の列も実装できる。
pub trait DigitalNet<const D: usize>: QuasiRng<D> {}
```

| 種類 | 変換 | 対象 | 保たれる性質 |
| --- | --- | --- | --- |
| `Shifted`（Cranley–Patterson） | y = x + u（`QmcWord::wrapping_add`） | 全列 | Weyl / R\_d の格子構造。Sobol では (t,m,s) ネット性が壊れるので非推奨 |
| `DigitalShifted`（XOR） | y = x ⊕ u | `DigitalNet` のみ | (t,m,s) ネット性。対象外の列は型で弾く |
| Owen スクランブル（将来） | ハッシュベースの入れ子置換 | `DigitalNet` のみ | ネット性 + 分散改善。v0.4 候補 |

- ラッパーはジェネリックな `Q::Word` に対して `wrapping_add` / `^` を呼ぶ。`urng_core::Word` にはこれらの演算がなく、u32 / u64 で別々に書く impl は重複扱い（E0119）になるため、`QmcWord` を介す
- `Rng` の Word と列の Word が異なる組み合わせも `QmcWord::from_rng` で扱う（32 bit の生成器で 64 bit の列をスクランブルできる）
- ラッパーは `inner` が `Seek` を実装していれば `Seek` を委譲する。変換は点ごとに独立なので、ジャンプ後も結果は一致する
- `DigitalShifted<Q>` 自身も `DigitalNet` を実装する（デジタルシフトしたネットはデジタルネットのまま）。将来の Owen スクランブルと重ねられる
- スクランブル後は「先頭点が必ず原点」という構造は消える。ただし座標がちょうど 0 になる確率は擬似乱数と同程度には残る（7.1）

## 7. 境界条件とエラー処理

`Result` もエラー型も持たない。不正な次元数はコンパイル時に弾き、実行時は panic しない（`overflow-checks` の有無によらない）。

| 状況 | 挙動 |
| --- | --- |
| D が対応範囲外（R\_d: 1〜16、Sobol: 1〜1024） | `new()` 内の `const { assert!(..) }` でコンパイルエラー |
| `Weyl::new` に偶数の α | 常に最下位ビットを立てて奇数に補正（alpha OR 1）。ドキュメントに明記 |
| 点を 2^B 個使い切った | 全列周期 2^B で先頭に戻る。Sobol は `trailing_ones` が B になったところで x = 0、index = 0 にリセットする（配列外参照を防ぐ） |
| `seek(n)` で n ≥ 2^B | n mod 2^B に丸める。`index()` も同じ値を返す |
| `next_point_in` で min > max | `urng_core::Word::randf` と同じく検査しない。ドキュメントに明記 |

- `const { assert!(..) }` は単相化後に評価されるため、`cargo check` では報告されず `cargo build` / `cargo test` で初めてエラーになる。範囲外の次元は `compile_fail` の doctest で確かめる

### 7.1 先頭点 x\_0 = 0

スクランブルなしの Sobol、van der Corput、s₀ = 0 の Weyl / R\_d は、最初の点が原点になる。`ln(x)` や逆累積分布関数で無限大になり得る。

- 既定では飛ばさない。Sobol は 2 の冪乗個の点で均等になるので、先頭を飛ばすとこの性質が崩れる（Owen 2020）
- 推奨する対処はスクランブル（第 6 章）。どうしても飛ばしたい場合は `seek(1)` を使う
- スクランブルしても 0.0 は出得る。`to_float` は u32 の下位 9 bit を捨てるので、2^-23 未満の値はすべて 0.0 になる。f32 で 2^20 点 × 64 次元を使えば、0.0 は期待値で約 8 回出る。`ln(x)` などに渡す場合は利用者側で 0 を避ける必要がある
- 将来候補: 開区間 (0,1) への写像（`(x >> 9) as f32 + 0.5` に 2^-23 を掛ける中点方式）を `QuasiRng` の既定メソッドとして足す。既定メソッドの追加なので破壊的変更にはならない
- この方針と理由をクレートドキュメントの冒頭付近に書く

### 7.2 型の性質

- 全型 `Clone + Debug + PartialEq + Eq`、列の型は `Default` も実装する（ラッパーはシフト量が要るので `Default` なし）。自動で `Send + Sync`。`Copy` は付けない。暗黙のコピーで同じ点列が 2 本作られるのを防ぐため
- `PartialEq` は状態の比較として実装し、テストで `seek` と逐次更新の一致を確かめるのに使う

## 8. テストと検証

「点が正しい」を既知値と層別で確定的に確かめ、「準乱数として機能している」を discrepancy と積分誤差で確かめる。テストは `tests/` の結合テストに置き、`std` を使ってよい。テスト内の乱数は固定シードの SplitMix64 を使い、結果を決定的にする。

| 種類 | 内容 | 対象 | 合格条件 |
| --- | --- | --- | --- |
| 既知値 | 先頭の点を参照値と比較 | 全列 | ビット単位で一致。Sobol32 は `scipy.stats.qmc.Sobol(scramble=False, bits=32)` の出力（2^32 倍すると整数として厳密）を `tools/gen_vectors.py` で埋め込む。D = 8 の先頭 256 点と、D = 1024 の先頭 256 点の次元ごとチェックサム（`SobolWide` は 21201 次元分を畳み込んだ 1 値）。Halton は分数で厳密に求めた ⌊φ\_b(n)·2^B⌋（scipy と照合）と 24 個のインデックスで一致 |
| 幅の整合 | `Sobol64` の上位 32 bit = `Sobol32` | Sobol | D = 1024 の先頭 4096 点で一致（漸化式は右シフトだけなので理論上一致する） |
| 定数表 | R\_d の α が奇数、f64 で計算した φ\_D^-j·2^B と相対誤差 2^-48 以内（f64 の丸めが j 回積もるため。誤った行の検出用）、5.2 の表と一致 | R\_d | 全 D で成立。80 桁精度での照合は `gen_rd.py` の自己テストで行う |
| seek 一致 | `seek(n)` 後の状態 = n 回逐次更新した状態 | 全列 + ラッパー | 上限 L 未満の n = 0, 1, 2^k−1, 2^k と乱数 32 個で一致。L = 2^16（`Sobol<1024>` は 2^10） |
| seek 連鎖 | `seek(n)` して 1 点進めた状態 = `seek(n+1)` の状態 | 全列 + ラッパー | 2^B 未満の乱数 n 100 個で一致（逐次更新では届かない大きな n の検証） |
| 周期 | `seek(2^B − 1)` から 1 点進めると点 0 に戻る | 全列 | 状態が `new()` と一致 |
| 1 次元層別 | 先頭 2^m 点が \[k/2^m, (k+1)/2^m) に 1 個ずつ | VdC、Sobol 各次元 | VdC と Sobol 先頭 64 次元は m = 1〜16、Sobol 全 1024 次元は m = 1〜12 |
| 2 次元層別 | 2^m 点のブロック（先頭から最大 4 個）が 2^a × 2^(m−a) の全格子に 1 個ずつ（t = 0） | Sobol の次元 1・2 | m = 1〜12、全 a |
| b 進格子 | 先頭 b^m 点が i/b^m（1 単位以内）に 1 個ずつ | Halton 全 32 次元 | b^m ≤ 4096 |
| スクランブル保存 | XOR シフト後も上の 2 次元層別が成立 | `DigitalShifted<Sobol>`（v0.2） | シード 10 種で成立 |
| discrepancy | L2-star discrepancy（Warnock の式、O(N²D)）を、一様乱数の期待値 √((2^-D − 3^-D)/N) と比べる | R\_d、Sobol、Halton | N = 1024。D = 2 では両列とも期待値の 1/4 以下。D = 4/8 は退行検知として下表の閾値以下 |
| 収束次数 | ∏(π/2)·sin(πx\_j)（積分値 1）を D = 2 で計算し、16 個のランダムシフト（Sobol は `DigitalShifted`、R\_d は `Shifted`）の RMSE を N = 2^8〜2^16 で求め、両対数の傾きを推定 | Sobol、R\_d | 傾き < −0.8（擬似乱数は約 −0.5） |
| 範囲外の次元 | `Sobol64::<0>::new()`、`Rd64::<17>::new()` など | Sobol、R\_d | `compile_fail` doctest |

**discrepancy の閾値**（一様乱数の期待値との比。実測は scipy と同じ点列で D = 2: Sobol 0.075 / R\_d 0.118 / Halton 0.141、D = 4: 0.201 / 0.762 / 0.275、D = 8: 0.672 / 0.824 / 1.092。実測に約 25% の余裕を足した値）

| D | Sobol | R\_d | Halton |
| --- | --- | --- | --- |
| 2 | 0.25 | 0.25 | 0.25 |
| 4 | 0.25 | 0.95 | 0.35 |
| 8 | 0.85 | 1.03 | 1.37 |

チェックサムは SplitMix64 の finalizer で `h = mix(h ^ x)` とする。XOR と回転だけのチェックサムは、デジタルネットの XOR 構造で 32 点周期に打ち消し合い、Sobol の全次元で 0 になる（方向数が誤っていても検出できない）。乗算型（FNV）も、先頭の点のように下位ビットがすべて 0 の値を区別できない。`gen_vectors.py` は、異なる列が異なるチェックサムになることを確かめてから書き出す。表の 1 次元を壊すと、テストがその次元を指して失敗することも確認した

N = 1024 では D が上がると L2-star discrepancy の差は小さくなる（D = 8 の Halton は乱数の期待値を上回る）。「1/4 以下」を全 D に課すと正しい実装でも落ちるので、D = 4/8 は性能の主張ではなく退行検知と位置づける。

**CI**

- `cargo test`（両クレート、u32 / u64 両方）。`urng_qmc` は両リポジトリを兄弟ディレクトリに checkout して実行する
- `cargo clippy --all-targets -- -D warnings`、`cargo fmt --check`、`cargo doc`（`-D warnings`）
- `cargo test --features sobol-21201`（全次元のテストに約 2 分）
- `cargo build --target thumbv7em-none-eabihf` で `no_std` ビルドを確認
- 生成物の鮮度: `tools/gen_*.py --check`（`gen_sobol.py` と `gen_vectors.py` は scipy を入れて実行）
- `cargo +1.86 check` で MSRV を確認
- `urng_qmc_core` の公開後は `cargo semver-checks` も回す（初回は比較元がないので対象外）
- ベンチマーク（criterion）: 1 点あたりの ns を D = 1/2/8/32 で計測し、README に載せる（v0.1 公開前までに追加）

## 9. ロードマップ

v0.1 で Weyl・R\_d・van der Corput・Sobol（1024 次元まで）を持つ urng\_qmc を、先に公開する urng\_qmc\_core 0.1 の上に載せて公開し、その後スクランブル（v0.2）、Halton と全次元 Sobol（v0.3）を足す。日付は決めず、各ゲートの条件を満たしたら次へ進む。

| 版 | 内容 | 次へ進む条件 |
| --- | --- | --- |
| core 0.1 | `QmcWord`、`QuasiRng`、`Seek`、`DigitalNet`、`AsRng` | v0.1 の全列と v0.2 のラッパーの試作がこのトレイトだけで書けること |
| v0.1 | Weyl、R\_d、VdC、Sobol（≤ 1024） | 第 8 章の v0.1 対象テストが全て通る |
| v0.2 | `Shifted`、`DigitalShifted` | スクランブル保存と seek 一致が通る |
| v0.3 | Halton、`SobolWide`（`sobol-21201`） | Halton の既知値と discrepancy、全次元のチェックサムが通る |
| v0.4（候補） | Owen スクランブル | 需要を見て決める |

v0.1 は最初の 3 列（Weyl・R2・Sobol 1 次元）から着手し、トレイト設計を最小コストで検証する。core は公開すると破壊的変更が高くつくので、v0.2 のラッパーを試作して `QmcWord` の過不足を確かめてから 0.1 を公開する。

**現状**: v0.1〜v0.3 の実装とテストは揃っている（クレートの版はどちらも 0.1.0 のまま）。v0.2 のラッパーは core に手を入れずに書けたので、core 0.1 の公開条件は満たしている。段階ごとに公開するか、まとめて 0.1.0 にするかは公開時に決める。

**公開前の残作業**

- Joe–Kuo の方向数のライセンス原文を配布元で確認し、`NOTICE` に転記する（レビュー時は配布元サイトに接続できなかったため、scipy 同梱の同一データから表を生成した）
- criterion ベンチマークを追加する

## 10. 参考文献

- S. Joe, F. Y. Kuo, *Constructing Sobol' sequences with better two-dimensional projections*, SIAM J. Sci. Comput. 30 (2008) — [論文](https://researchcommons.waikato.ac.nz/handle/10289/967)、[方向数ファイルと sobol.cc](https://web.maths.unsw.edu.au/~fkuo/sobol/)
- M. Roberts, *The Unreasonable Effectiveness of Quasirandom Sequences* (2018) — [記事](https://extremelearning.com.au/unreasonable-effectiveness-of-quasirandom-sequences/)（R\_d の定義）
- A. B. Owen, *On dropping the first Sobol' point* (2020) — [arXiv:2008.08051](https://arxiv.org/abs/2008.08051)（第 7.1 節の根拠）
- B. Burley, *Practical Hash-based Owen Scrambling*, JCGT (2020)（v0.4 の Owen スクランブル候補）
- R. Cranley, T. N. L. Patterson, *Randomization of number theoretic methods for multiple integration*, SIAM J. Numer. Anal. 13 (1976)（第 6 章 `Shifted`）
- T. T. Warnock, *Computational investigations of low-discrepancy point sets* (1972)（L2-star discrepancy の計算式）
- 先行実装: [sobol\_burley](https://docs.rs/sobol_burley)（Rust、Owen スクランブル付き Sobol）、[scipy.stats.qmc.Sobol](https://docs.scipy.org/doc/scipy-1.10.0/reference/generated/scipy.stats.qmc.Sobol.html)（既知値の照合先）

リンクは検索結果で所在を確認したもので、本文の精読はしていない。書誌情報は実装前に原論文で確認すること。

## 11. 改訂履歴（レビュー反映）

主張はスクラッチのプロトタイプ（rustc 1.97、urng\_core 0.1.1、scipy 1.18、mpmath 80 桁）で確認した。

| # | 変更 | 理由（確認内容） |
| --- | --- | --- |
| 0 | 第 2 章の前提を訂正し、urng\_core の修正を公開前の必須作業にした（0.1.2 で修正済み、依存を更新） | urng\_core 0.1.1 の `u64::to_float` は最上位ビットを捨てる（`(1u64 << 63).to_float() == 0.0` を確認）。旧版の「上位 52 bit を使う」は誤り |
| 1 | `QmcWord` を core に追加し、`QuasiRng::Word: QmcWord` にした | `Word` には `wrapping_add` も `^` もなく、ジェネリックな `Shifted` の `next_point` が書けない（E0599）。u32 / u64 別の impl は重複扱い（E0119）。`QmcWord` を介す形はコンパイルを確認済み |
| 2 | Rng の Word と列の Word が異なる組み合わせを「非対応」から「対応」に変更 | `QmcWord::from_rng` で変換を書ける。sealed は障害にならない |
| 3 | 1 次元列の `Rng` 直接実装をやめ、`AsRng` アダプタにした。目的 1 の文言も修正 | 連続値を組にすると π の推定が 3.29 になる（4.4）。同時に `Q::Word` の曖昧さ（E0221）も解消 |
| 4 | Sobol の V を `&'static [[W; D]; B]` にした | 値の const だと debug で 512 KB のフレームとループ内 `memcpy`（D = 1024） |
| 5 | Sobol の次元上限を feature なしで 1024 に、`sobol-1024` feature を削除 | パラメータ表は const 評価専用でバイナリに入らない。v0.1 の Cargo.toml に v0.3 の feature があった不整合も解消 |
| 6 | `sobol-21201` に `long_running_const_eval` の allow と `fill_point` を追加 | 21201 次元の const 評価はリントでエラーになる。戻り値 170 KB |
| 6' | 6 を撤回し、`SobolWide32/64` を feature で足す形に変更（5.4） | 実装してみると `#[allow]` は効かなかった。リントは利用者側のクレートで判定される。チャンク分割と遅延評価で、利用者が何も書かずに全次元を使える |
| 7 | Halton の 64 bit 計算を u128 累積に変更 | 「最大 64 ulp」は誤り。上限は (b−1)·K で基底 131 は 1300 ulp、実測 559 ulp |
| 8 | Weyl の状態に `seed` を追加 | `seek` が s₀ を使うのに状態になかった |
| 9 | `Weyl::new` の偶数 α を常に補正に変更 | debug と release で挙動が違うと、release の経路をテストで通せない |
| 10 | インデックスの加算をラップ、範囲を \[0, 2^B) に固定 | `seek(u64::MAX)` 後の `index += 1` が `overflow-checks` 有効時に panic する |
| 11 | 「スクランブルで先頭点問題も解消」を修正 | f32 化では 2^-23 未満が 0.0 になり、0.0 は擬似乱数並みに出る |
| 12 | seek テストを逐次比較（小さい n）と連鎖比較（大きい n）に分割、周期テストを追加 | 乱数の n を逐次更新で再現するのは計算量的に不可能 |
| 13 | discrepancy の比較先を閉形式の期待値にし、閾値を D ごとに設定 | 「1/4 以下」は D = 4 の R\_d（0.76）、D = 8 の全列（0.67〜1.09）で正しい実装でも落ちる |
| 14 | 収束テストの被積分関数と手順を明記（ランダムシフトの RMSE） | 単一の列の誤差は N ごとに大きく揺れ、傾きが −0.11 になる例があった。不連続関数（円の指示関数）は理論的に約 N^-3/4 |
| 15 | Sobol の期待値を sobol.cc から scipy（bits=32）に変更、`Sobol64 >> 32 == Sobol32` を追加 | sobol.cc は 32 bit 精度の double 出力で 64 bit 版を照合できない。配布元にも接続できなかった |
| 16 | 2 リポジトリ構成に変更、依存は `path` + `version` | 実際のリポジトリが `cet-t/urng_qmc_core` と `cet-t/urng_qmc` に分かれている |
| 17 | core の 1.0 は `urng_core` 1.0 が前提、と明記 | `Word` / `Rng` が core の公開 API に現れる |
| 18 | 孤児ルールの記述を修正 | core での `Rng` 包括実装は「衝突する」のではなく E0210 で書けない |
| 19 | `const { assert! }` が `cargo check` で出ないことを明記し、`compile_fail` doctest を追加 | `cargo check` は通り、`cargo build` で E0080 になることを確認 |
| 20 | `DigitalShifted` も `DigitalNet` を実装、全型に `Default` | スクランブルの重ね掛け、clippy の `new_without_default` |
| 21 | テストを `src/tests/` から `tests/` に移動 | `no_std` クレートの `src` 内では `std` を使うのに `extern crate std` が要る |
| 22 | Joe–Kuo 方向数の著作権表示を `NOTICE` に置く | 方向数は BSD 系の著作権表示付きで配布されていると記憶している（配布元に接続できず未確認）。MIT / Apache のクレートに埋め込むなら表示が要る |
| 23 | チェックサムを SplitMix64 の finalizer に変更（8 章） | 初版の XOR/回転型は Sobol で全次元 0 になり、1024 次元の照合テストが何も検査していなかった |
