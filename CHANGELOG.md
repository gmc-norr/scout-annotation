# Changelog

## [0.5.0](https://github.com/gmc-norr/scout-annotation/compare/v0.4.0...v0.5.0) (2023-05-02)


### Features

* add AO and AF to sample FORMAT if missing ([813876c](https://github.com/gmc-norr/scout-annotation/commit/813876c33d65007b5de659c303a51a36c940ff32))
* add bam file to load config ([bebc322](https://github.com/gmc-norr/scout-annotation/commit/bebc322bd0d2a0e987419542e51bf4ede1c67c5c))
* add notemp flag to cli ([f9d763d](https://github.com/gmc-norr/scout-annotation/commit/f9d763d2ecc158ee3d595699bd0fad45db3612f9))
* add variant normalisation ([a94b623](https://github.com/gmc-norr/scout-annotation/commit/a94b6231c60e220838f88775d9b4411841b6fbf3))
* add VCF filtering ([87646d9](https://github.com/gmc-norr/scout-annotation/commit/87646d9fa311a9d878392033387280a493b6e649))
* possibility to skip filtering steps independently ([38f1138](https://github.com/gmc-norr/scout-annotation/commit/38f113821abcbdfb6c65a4cd6b589c97faa1f2e1))


### Bug Fixes

* add sample name to log output files ([8f0f692](https://github.com/gmc-norr/scout-annotation/commit/8f0f692eb63063c05d8ad572665254abe190eb16))
* allow missing panel and filtering columns ([4ddcc70](https://github.com/gmc-norr/scout-annotation/commit/4ddcc70aff9c7b64539a3fb7fbb07ddfa19742d2))
* AO filtering failed for missing values ([105603e](https://github.com/gmc-norr/scout-annotation/commit/105603eb36daaca463ce3c96d90952f072468850))
* better checking of sample parameters ([daa1885](https://github.com/gmc-norr/scout-annotation/commit/daa188514433c44b385d3db1af49a8e9866ee9b1))
* check that the bam column exists ([760f4cb](https://github.com/gmc-norr/scout-annotation/commit/760f4cb29028031f326a182914fe3c9b146e37a1))
* correct header order for vembrane ([a3c33df](https://github.com/gmc-norr/scout-annotation/commit/a3c33df7e60f14e2b3a020da9aab2473c006fbf1))
* don't run genmod if there are no variants ([e27fb78](https://github.com/gmc-norr/scout-annotation/commit/e27fb78f8eb0780ac3f3d8be10a89501db3990c6))
* keep uncompressed vcf ([f4691ba](https://github.com/gmc-norr/scout-annotation/commit/f4691ba09e5273cc51d3484e6870f52c4f46dbb1))
* missing VAF for variants with faulty AO ([0809f57](https://github.com/gmc-norr/scout-annotation/commit/0809f570ecc2302a094f416f4b60581207b0693a))
* rerun rules for incomplete files ([500c76d](https://github.com/gmc-norr/scout-annotation/commit/500c76dbd4e9b7ec0870cfde90e2bc48bfe52e7a))
* split on entire file name in batch mode ([fb5a58c](https://github.com/gmc-norr/scout-annotation/commit/fb5a58c09d20d0a964ad9bd3d2b168c938a96e40))
* update package requirements ([6fa7e4b](https://github.com/gmc-norr/scout-annotation/commit/6fa7e4b35439f589b857f0932111f1affdd35edc))

## [0.4.0](https://github.com/gmc-norr/scout-annotation/compare/v0.3.0...v0.4.0) (2023-03-15)


### Features

* add batch command to cli ([5f70d29](https://github.com/gmc-norr/scout-annotation/commit/5f70d29e5b714187385951ed0fa0571b41f7958a))
* use repo dir if local dir is nonexistent ([ccb01b5](https://github.com/gmc-norr/scout-annotation/commit/ccb01b5df56bc233fb646155dc676811be438a29))


### Bug Fixes

* better error for missing panel ([9641d00](https://github.com/gmc-norr/scout-annotation/commit/9641d00b44d6d293002eef4f647d46ba37ceef82))
* get version from command line ([d57a285](https://github.com/gmc-norr/scout-annotation/commit/d57a28522e0d563a63bd26d83ed6128222bb3a54))
* relative path not correctly resolved ([f4f5d2a](https://github.com/gmc-norr/scout-annotation/commit/f4f5d2ac7c58f63a2d4601e8d4cd0f5497ce42ab))

## [0.3.0](https://github.com/gmc-norr/scout-annotation/compare/v0.2.3...v0.3.0) (2023-03-13)


### Features

* add command line interface ([ef453b4](https://github.com/gmc-norr/scout-annotation/commit/ef453b4ac6d1b141d65bbbccdc20089f1862d92d))
* add panel filtering ([d369853](https://github.com/gmc-norr/scout-annotation/commit/d36985385af90c5cfe1fe8b7bcdfa82825dbee80))
* add panels to load config ([6a94e14](https://github.com/gmc-norr/scout-annotation/commit/6a94e1476b921bbcf2c01eba21d3f80a3a65b699))


### Bug Fixes

* better cli path handling ([68ecdd6](https://github.com/gmc-norr/scout-annotation/commit/68ecdd64e7ab595002759e72975784525177d846))
* better sample file names ([51fa1f8](https://github.com/gmc-norr/scout-annotation/commit/51fa1f859d66a532f873e292fc72d551f57eb6d6))
* corrected panel metadata ([44fc063](https://github.com/gmc-norr/scout-annotation/commit/44fc0633aa9df1123b0cdd6b01e7d6a4978bd74f))
* don't filter on an empty set ([b5dbe27](https://github.com/gmc-norr/scout-annotation/commit/b5dbe27506baefa7daba0746c0e42a87b9b7fae8))
* snakemake linting ([edeaab6](https://github.com/gmc-norr/scout-annotation/commit/edeaab696cdef2336f3d9db1625547fbd826139e))
* update paths to config files ([4c4fe0e](https://github.com/gmc-norr/scout-annotation/commit/4c4fe0efad697aa6952555be3b5498f9af1ee98e))
