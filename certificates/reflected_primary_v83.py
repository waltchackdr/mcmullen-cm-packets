#!/usr/bin/env python3
"""Full-only V83 portable replay frontend; generated only from reviewed data.

POINT_PASS means one fresh complete target/label replay matched its frozen
expectations. It does not adopt executions, register a theorem or finish a
release. No preparation/source-only mode exists. Failures retain all output.
"""
from __future__ import annotations
import argparse
import ast
from collections import Counter
import hashlib
import importlib.metadata
import importlib.util
import json
import math
import os
from pathlib import Path
import platform
import re
import resource
import stat
import sys
import tempfile
import time

ROOT = Path(__file__).resolve().parents[1]
BACKEND = 'primary'
KERNEL = 'certificates/reflected_primary_kernel_v83.py'
SOURCE_PINS = {'certificates/d7_ambient_points_v45.py': '34cb35db91b953162b595fb61e031f551e6df0a0d7b5192dc15a6c00c39274a2', 'certificates/d7_r83_mask14_v46.py': '70e9c7114adea2b319c3e430046ae3e917717142e7a78e4d43a043950f1147f3', 'certificates/d7_r83_masks14_30_v47.py': 'e78ba5e793134f28d8113d8129809421a5a7754c439afcd64ba0385543799941', 'certificates/d7_r83_principal_expected_v55.py': '0fd36a67b004e6a250188aa8731a207456b5e49faa955005d6b12cd36665bea7', 'certificates/d7_r83_principal_v55.py': 'a937b69037b5d53f40a600aefb177bf89766d3529b22f73529b2b861ec9238cc', 'certificates/full_packet_v62/ambient_points.json': 'c0e84da6c266c087cbce65a6e8aca8e3f30350ce56b5a1089c4ae4978a381f36', 'certificates/full_packet_v62/module_templates.json': '09ccff217511a67e946567135f9ffb74da020843621d0895ecce8746d25c245d', 'certificates/reflected_packet_v83/expected.json': '0ca445d277730d443ea5f37ddfde362836ed591edca177ed20d8dfcc45e4eee0', 'certificates/reflected_packet_v83/plans/original_q_v1/om-pic0-00.json': 'c6421803dfe558b66b18c01526b33c1c734d9718beca5e61c8620e8a35f3a83e', 'certificates/reflected_packet_v83/plans/original_q_v1/om-pic0-07.json': '601e882858830d2cb70f1b5a3102b1acfe1f7b6ec04095a08f25104f318400f2', 'certificates/reflected_packet_v83/plans/original_q_v1/om-pic0-15.json': 'e685879b3ec800514ba0fc8043133ee0ceef83965fa16a2b16238705acd0d423', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/om-pic0-00.json': '693f4669ef550d4a15d8b5c96849a16294e94693fbfc36c60821505136e73147', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/om-pic0-01.json': '836dd0ad28d48b9128a5bebafcd1f03de8de2e0a96b6180238c52109c2cf7100', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/om-pic0-02.json': '11c7c8ccd8b84947bb332e2df01612ab84c36a9adcfe01a4ee2fd9440d47e1f2', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/om-pic0-03.json': '06dc85ae6fd81ad8e077736b9a99f58ddcf3936a75e4f47cf7cab42b65daa172', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/om-pic0-04.json': 'af4ad29ed59ae9afdd5d8f539e95ca6ee6548be3992e683f0fd6336c96bbabe1', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/om-pic0-05.json': '158b0e626475ee3f95eb0dae3c4021ea7ad476136c60960f0a89bad5fd4d186d', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/om-pic0-06.json': '1c2c627b38f0ceae8a2eb9e33801763b019925de46618d057fb15dfec9565966', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/om-pic0-07.json': '4f80641a798f928c574afb3263a0fb4a6f1f300c4a7225b24e79e523005f6a86', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/om-pic0-08.json': 'f299ea413a746067b3060cde9ab658b038de8130fa39a8d3d16263d9a4fb0533', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/om-pic0-09.json': '5f2dfb998b63ac6187df6bb56571fed47926f55536740f23f5ddc8bae62a83b6', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/om-pic0-10.json': '92dbb6b19f24ab6ba0e9ca9a4a8a7aeb99d8c6578c254ea83028d37d27cfab00', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/om-pic0-11.json': 'bca4438617836dc87c80c71d470a6ecab57af1017a186eeaecc17eaccbea1cd2', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/om-pic0-12.json': '68aece371f51bebed19941c670efee1d9a107699176c3c4be2a8015e62c80124', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/om-pic0-13.json': '5eaa79d08adcf5f0c813aee3cc366196c567a4b9ff6e9412b529cd92aae3aa91', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/om-pic0-14.json': '8b41e2d37fc84c63f5a8c9f35ba0803b2b96ce9325cc51e2ac60ad6815452d27', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/om-pic0-15.json': '4baa75cbe0065eacdec6fbe43e55687b15de0ff9f3bd86b2da32704b410f3aad', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r127-pic0-00.json': 'f51b9f93f71836e2f8ac622dc46da2b81f169a93f2e07932f9d6981b7d9f686b', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r127-pic0-01.json': 'e5f5cc9f5c7cf9e9d33c05e961c44099445d6890c705c6f71c2821fccb7a1c73', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r127-pic0-02.json': '2ffaa48b050a4d51220038da826fe962e5ee0671b29e341bcd763380ae4216ad', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r127-pic0-03.json': 'eb7451fe57e16c7ade3b8a4a285df75650a41ea57b7c6e7fb114bd8f7d7f81d4', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r127-pic0-04.json': 'fdde636ddde018df2a57e63bcd04827e4866c42da74a8a9fc06794def9413d32', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r127-pic0-05.json': '6c4551316f7fcbbbbffdeb14cce83894d240f1357c7bfad8b0ed6b4d3de63c9b', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r127-pic0-06.json': '4a16efcde635d782b5788c12777286f974189b90f76861bd1886b42a36b5cee8', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r127-pic0-07.json': '4186175eb5875fd92d1b38e03ee51de4921c491c0ed1c618176dcd7e0788c06c', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r127-pic0-08.json': 'd16525db3ff39e7a814c06a2a91651ab12cee9c5330bf93ec260de444d63d9a1', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r127-pic0-09.json': 'b9d656666cedd9539ae6fd2a05a9c66d91f4b75ba342ffef5fc984254496294a', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r127-pic0-10.json': 'b045881974d367c4b81ccd402cb7c8764b795b3ab87fa2a49b8678dcfb580d2f', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r127-pic0-11.json': '90973eba2c79fe578a75c313ccb22d3b724337e7397cf595bfbd63b0b3a9b8d0', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r127-pic0-12.json': '6c032e7205b92ee1d88354e31f0a4e232ba8f559f8ebbf7447f13d19bc953481', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r127-pic0-13.json': '9f34f17f74a84e220c979b1461afb93cd73337f7ba0c2eae4e44d3ee60e76268', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r127-pic0-14.json': '47dc789c51d6b1cee01aee7e5ddad77c3b07243087d0f7b4903f2aa44574833d', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r127-pic0-15.json': 'b958f62a088a2f683b34c4c8be1b456c6be7c20d3ab3c375d53164075b491c8f', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r127-pic0-32.json': 'cc15ba4d57666e7cef3a6c17a57fa201779b54d47caa666f73583d26a1a9d101', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r127-pic0-33.json': '918448e64013d09e16a77511df4282ef5f1a1691ec03f45907a794e3db701b36', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r127-pic0-34.json': '7604e85312f2f0f78b8992c42c00be5463cb6489f527c76e47ed5753c44ac90b', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r127-pic0-35.json': '71f4528c0a178971d6ea3606c6630b658dfbcd17b3e833b81b9688fbf6d2782f', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r127-pic0-36.json': '5bd36b24d3810b27cb43cdb8cd8192e5e995152f5025cff3d0cf0d4a2169cd83', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r127-pic0-37.json': '3c0f0ab08e6db6b67f5c8952d118dacb3b7a80534f01dca4ef2b8bf2f27f45a2', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r127-pic0-38.json': '4f928a30219032918da872598168814ce7b398e6264aebe7dbf506050f8cfe34', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r127-pic0-39.json': '18869b9584c3bf9ffaff03386d7bad24e95b59736d95a6ee3e16b62b31a043b5', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r127-pic0-40.json': '41f34364a7343fd9add10366ceebacac05c39314d62f5078fc49182749a2ccfc', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r127-pic0-41.json': '7f6c5ddf51622e865ca5e9fe52b3c9caa74ab338fbbefd0cc3b9d51aab0d52b1', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r127-pic0-42.json': '69228266a77125f9560d17d9a232784aeb103fb433315ce2515bef8abdb01846', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r127-pic0-43.json': '3680d95f4c09f96db07ff86a95ab417f2e856944375f3bc7c69ba49c6818d768', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r127-pic0-44.json': 'd16c15bafc4cfc197f710c705fc296f9c7ec4d80a833b1118cdb6d0705f5619b', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r127-pic0-45.json': '8b2c2aa0faca176f6eeeee89068648e673a655adfc1c142c40d9196dea8e7f20', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r127-pic0-46.json': 'a9d96e0f31294c74cbea51e79e7e6de692a9bb76dd0dfd31d2afeb65a38ac1cb', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r127-pic0-47.json': '9f05168a389931ba3eee59686b2313fdd05f22bfc02ab229d09d8d49fca842de', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r127-pic1-00.json': 'b93ce79bbb184436cdb1da6d2c4e579e8b9be7bf5462663589fd358403002278', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r127-pic1-01.json': 'ee9d3bcab87638ea4481749284a0497adb65d82b782170498b65459dcab4ff00', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r127-pic1-02.json': 'd6ac92c3715724b977fa3065c3e95d7e4fa2885a9c48c7b41f3b0722bf77fd06', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r127-pic1-03.json': 'a01f2b2f07e2766d3f735923b353020f91c6869de5beae41122bf74643225230', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r127-pic1-04.json': 'f11654e08419f32e15ad539d4cfe4259b6cd0fe56fcbb6a0c70dc9823ded9dd4', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r127-pic1-05.json': '81c686be900debe1475154539455681771ec4a1e7efbbf2ee9667372483aab28', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r127-pic1-06.json': '0775a6619121da7528171e30ba95cdbb56eda298f76eddd90e652b9ca3bb2217', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r127-pic1-07.json': 'ba566a8278fe24bbb09700ff90c6d882b768158b4eb7ad7b3b06ec1462c4fe24', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r127-pic1-08.json': '0607cc05cd05906e700dac2659daa623f8620513278498b20c92d5ee52201594', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r127-pic1-09.json': '75d172eff1562bc7b7241425f6ed2d1ab9f714da22663626d3ce769a75c6a760', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r127-pic1-10.json': '70289d8778f54a44b893957b8b1cd235fdb57b6a0a92fe5dd70e1c4e1172d15b', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r127-pic1-11.json': '4fd1031516b4f2a97485590ddcdd4cfd0b662cf4dbde97ac87b282cf44c727f8', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r127-pic1-12.json': '99a362965a4ecb89e9fdcd8d082a43151be411b6200914fe558f940a10135f52', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r127-pic1-13.json': 'dfb88414b08c50822b903204389ea5f5818e6ff9aa7e6a151fb41425b61f5d35', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r127-pic1-14.json': 'b24ceb85ed056753547cc245c8ca18b4757a7e4dbe57449333cc6f7cd6cc281f', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r127-pic1-15.json': 'a298ed23850fd2fba32287e9c5545fd984a33c567cebe55b5b8c9b11e4619be7', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r127-pic1-16.json': '7faad2e6bd33c9cdfbceb8087b2d35cf0f19cb04d940d01f2e53f9f81ecee699', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r127-pic1-17.json': 'ae0e6d51e45da79b32372ff4c9d9b84e93c096b2d5fe90d7fab290a6534b06b4', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r127-pic1-18.json': '04527b502c2be9a87de5ab152000bacd8b3777e39a5aa08dc4b772b3dea3b283', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r127-pic1-19.json': '06c814cccc7f212843c5c4b9a1908bccabd9206c3bfd4a89d563539b6e06fbe5', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r127-pic1-20.json': '5e235007d38b13f209c10d7eac94fa7cbcaff8d1ac26045943c4b2dfa53b021e', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r127-pic1-21.json': '56d32a1d88eaac7d3e483de650753532386ea13ec3ab0190bfc26322693221c5', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r127-pic1-22.json': 'd3f1ffef3a3d8c2f0a529289ee8156408088c114745dcd297358067a92beff8e', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r127-pic1-23.json': '7bdb2ee788a280c15b1aa5a32c7f734336d3d8e3030c2c184cccdc257ddb40f2', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r127-pic1-24.json': 'ed3153523cfa132eeee954ab3e9c17585d408fd582fc7b9af3116a32822a8fe1', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r127-pic1-25.json': '6bb5bd07c743fd16074aaae8b29a24a9c80196dbd9d4a09959a4e6e7a0843c8a', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r127-pic1-26.json': 'ab15feb0625a3e74ed10eeaccba52eb6e25b7b2837015c177442714525cb58fd', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r127-pic1-27.json': 'ca34185ea857f85611bd0dccde1eff564d44c9100832b56310e37af705c556ba', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r127-pic1-28.json': 'e83fc5df4a762ac4002bc1b066fb71c6880a8dbb5cb379b7a08684213f6df3a3', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r127-pic1-29.json': '05f8449a62480c4cd121660c5f5d7841dd260d99227b9c56c027c6339b3ac77d', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r127-pic1-30.json': 'e097ad51fbb6eba82afbe1965b9ddc2df5b445c315619dc41e0ecdffdfce8f53', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r127-pic1-31.json': '48b144663398cc489c088e86605081f060b2132060147b3c5238e67067d4ae53', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r127-pic1-32.json': '1d245236688b296eeb228d014d76bbc9a489cfc089c15bec99f6d8cda8d83058', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r127-pic1-33.json': '64b765d142f715ed529f77b79561c0a23e4c3cbe7cd5b6fe5a3502957302328b', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r127-pic1-34.json': '2e23e4b6c2ee8b83cd3e821f0cfd54ab9dffca60d668fdb3542dbba12eda2352', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r127-pic1-35.json': 'e2694c2c4237d0bef0399b4e79c18f6f4bf4ceb932dac72660bebe42c51a2de6', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r127-pic1-36.json': '951c454d21e387603b45574b125519f584b04f6692ad994d7dfa1927e5778655', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r127-pic1-37.json': 'c7883485b82794e7fe31d2f952b4acfd85c3e7094154ec084554630723d9e123', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r127-pic1-38.json': 'b1b805dca18b74ecd33f27b5abbd06fd5c01b74cf42ff095b82c647782a694a5', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r127-pic1-39.json': 'be6c3dd34628054522b6b8b74d5d2ca1852f62886b367f8d88764ac58781ff45', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r127-pic1-40.json': 'abb3888df5511902ce2eb8c86fde4c5004485142b9ff9c2c5377988b91ebeafa', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r127-pic1-41.json': 'b346905e2dbdd726fd2a60ec1f3fc72d13b9bc043563f991465e977051d74bfc', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r127-pic1-42.json': '3d3f1de32068ba6a3851e0cf56dc05e8382ceb30e26ca31b766f3037794a5151', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r127-pic1-43.json': '1056477e4cc606a52759736e4d68dd203e9df20f70a149d8f114256af34814cf', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r127-pic1-44.json': '8fd6013a83916443b241641081525dc9cfc2a3791c4527c418444ba292921099', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r127-pic1-45.json': 'ac9d4cbf8e8303895252795dee61b7bf0ac346bd3b99dac2aee905f9a02e21b0', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r127-pic1-46.json': '28df0ea2919b8e273792452acba34e39ba1710efd2ce581b95524b42beb2d1e9', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r127-pic1-47.json': '55a2be2d08df3b84aedc0b5455ebb99c470b891b56aaea662b5094d6522ed77a', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r127-pic1-48.json': '853e6a98ed6616ef56b7236a0875d190b4f404f221bd46c28f2ec9600be6535e', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r127-pic1-49.json': 'c9bcd67c2ac76a1fe357f21d59ae0e471c2cfaad3d6c706eb6a8f5129deb0311', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r127-pic1-50.json': 'b17b712ea80a66d4255aed65eecc599acb25d5d9871517411960b8275142dd5f', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r127-pic1-51.json': '2e95a49a99146738241c0422ecfa86c98c9bb33708c62c26f90689917c0edab8', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r127-pic1-52.json': '05b26e16a4c0d8444ce77cbdb7290478208b14f260cb1771435071121fd3a102', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r127-pic1-53.json': '300bca7be4f244618c4f2be1d88e7853f8277527864a96fd23bbc56e5ab408c4', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r127-pic1-54.json': '7a735e953d62b4ee0db976a0dc07f76c8105fec4ea338289fbdb7c6230c8a50e', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r127-pic1-55.json': '481106ccc373d85be3b88b0c9aa4c29023db16265d08aa56a7294233da560f7f', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r127-pic1-56.json': 'a16fb6302f2eaf1ca179a18adc740a144be867172d2b9381a90021a8433ee029', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r127-pic1-57.json': '0616e0fe1752a2564a269563ea28620347c6f4939def116cb4ed47492da83024', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r127-pic1-58.json': 'a0c8baeca26fc75b812de031a6ebb116160559391080107fea851994baff62a5', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r127-pic1-59.json': 'f2737a0bf932a6be16b2cd354b949d288b8dc3cda755f197092a00ebbc007ad9', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r127-pic1-60.json': '4a27a274c687806b7919ee59b7aaf6f19f8174f11a58521088cf9bd094fd947a', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r127-pic1-61.json': '7edc353d90372581c0541bd787f7c7bbb4ea41eeccb3af8bb37fe48c842bf64d', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r127-pic1-62.json': 'ac7d3c34c1adb97eab9d59aeab61c4a76f05d322ab85ac4f26c17ba890ede2b4', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r127-pic1-63.json': '2d3e9779de5be3f829636fb4ffe5f7a1fe93e156d83b53938bdc467a4b7d277a', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r41_12-pic0-00.json': 'e1850c03501e0a1e062a588b5d0cd84df0bf7b6491493a0dc92b174002bd7045', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r41_12-pic0-01.json': '5a6646497e961a5c7e30d38d6125a6c8dadbe22dd71dcc102a00526e70659d23', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r41_12-pic0-02.json': '5c80bdd5dd2afd4f8d7bbce5b76d18355469d7e5f58a72bc6432cb813f467454', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r41_12-pic0-03.json': 'c0e11eda3b4a0e654ede1bd5c2020d0bc48ad74682a55f1b4b352277ec88c5ec', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r41_12-pic0-04.json': '107b4ab0e64ba065951478ebbb28c54ce178b3eed4f4382df12dea48b3425614', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r41_12-pic0-05.json': '6a0ddabf36a8d8c952b73b426eb66fedd6962c78b108675c8606af185798d280', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r41_12-pic0-06.json': '12b97a29ec03198c4c23cfaf6c9fd5f74a2d114e93767097c0a0bf395a3274f7', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r41_12-pic0-07.json': '79602b73773b4d0d2b8d7ca799bd9dce1441fa4777b37a75cdd82d8cde29ae3f', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r41_12-pic0-08.json': '281a9fe85ef02edfe0826956b37db2b4c9f9427342dececc2293d08bde7d766f', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r41_12-pic0-09.json': '32f07043fce59420a84eee3184263ec9872febc7eaba1185a1487b348fb6713a', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r41_12-pic0-10.json': '523a48864ed218a384aa39b63a7ae7b92936e318cb948be8a52b6737017b94c1', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r41_12-pic0-11.json': 'b9732d38c0b4f1e9e3344330d49af11ceb9f7e85c61d6ffe2aa45b228e61ba60', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r41_12-pic0-12.json': '773f0447cae07ceb7cc876e77673f3b6e375d58f0f09fcaae77ef877f29c13b6', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r41_12-pic0-13.json': '40bca25786f72d5f85899085647e244445db0bd63122a82752177c213e2087ae', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r41_12-pic0-14.json': 'c0d4451bd08e18512aae67e8ed2d91cbf5dfb09d1643b74f86c030b3bdcefa84', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r41_12-pic0-15.json': '64e8dbd6e45eb8b36652feac6d2df5d4a05fc6aa2c27056aa845740d7586129f', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r41_12-pic0-32.json': '1dc6f96aabc571da247dcb4aeb10f5090a2c7fdebd97c92d48a693dc5d0dc38d', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r41_12-pic0-33.json': '4d5a5ac889c1e5d6b8d135ba2721d45f46a571303310129b52e8f099e24cb4e3', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r41_12-pic0-34.json': 'd01dc0b76cc429a647ab5f95de323cc4a126b973d6e23ebad27f164d88d52ed0', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r41_12-pic0-35.json': 'b257b2da7937ccb08a53ebf8a2496560a597cc4f609ddbf768b53c25252c55f8', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r41_12-pic0-36.json': '2ae097bf447909cc50cf194252827272e251f16067453a9928e29d5247053edc', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r41_12-pic0-37.json': '163764f994746022de41dafd2726440bbcf113a584c425793fafbd44607fc566', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r41_12-pic0-38.json': 'c74e2bfb713ea1c0dc0e09976f06ac579fe449a0ee5e779b9840068f7be4ca30', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r41_12-pic0-39.json': '76a78dbedf83aa36458e63cd9fdc58ae6d6efd6b5df6c2e652f27a737673a84b', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r41_12-pic0-40.json': '3b19575c1e373e7163dde4baf9fb699918b95e98b8cc92fec3f2ef9fa33ab25f', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r41_12-pic0-41.json': '0c9001d8405f8d2b7eef227ff545b9e48ffc28da2510fe96af22bfc2dbcc5737', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r41_12-pic0-42.json': 'd54872108425b725783ab2baff96a322978d3c2087f2270f092d40945be6fa42', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r41_12-pic0-43.json': '879a38c25be7447934f469e3838343c37bc158a15827bdf2168beec20fc7b0e3', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r41_12-pic0-44.json': '2c5814324bf8e111cbf0b692d144bffd9f3b42873d77e4eef18fab0f3dbbd708', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r41_12-pic0-45.json': 'e3e14ba392c20f9d8d3a7fccfb13a1ee6ae2d721b420ab4f3b2733eaacaf9127', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r41_12-pic0-46.json': 'b33e1fcc0897c0825a0bcf4f2eee10352760070440a6b770161a3f1e95b8bb1d', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r41_12-pic0-47.json': '47d010585cfa4bb6ea69078aecffd55f9e5e0b82de811c2b26bd45e6ffb90911', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r41_22-pic0-00.json': '23d48ed3db898fa17079ba071f1eb1932ed223a91cccbdb8f1c77bb2bc6f80a6', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r41_22-pic0-01.json': '27c1fc855f7afe90d506512cb8d0f62f707d67ef1b25d9ced4a37f45cf8b7ef3', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r41_22-pic0-02.json': '100a1442afc275785f1a5bafff078aaf53ecee2ce7c89ab63ed123c9015f5ce1', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r41_22-pic0-03.json': 'b7c7f08c9aaf0a82c74e84b401e57235c6a678ac2d53e65d10887c6d0d6065e2', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r41_22-pic0-04.json': 'd900f103b8225b3252162c79c93f4b12655649932999a947665534d0c259057a', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r41_22-pic0-05.json': '52c7d3cccc6ef36e5b85dd3425385e85efa7c4b4ca06f16e29dda96fbedc75c1', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r41_22-pic0-06.json': '213a8370fb811d74ccb879ba23193dd50f817c8ce66a0a2671b257560f4c90ca', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r41_22-pic0-07.json': '24f046acd2cd89132207a5618f881e2d26f848ff72e9ef89aa506bab1658368c', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r41_22-pic0-08.json': '6dc47aab0aad9174c4272c23d5af4b204b426fbd0c0d90eff6b7fc5dbc4e96ad', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r41_22-pic0-09.json': '616f13f57805924fe6c8dc27fbdb10d1f711fd5d571b58a005a22323edcceecd', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r41_22-pic0-10.json': 'f816b508d64f23863c781d1890d421370c19dee83403255f4c47417ca46691ef', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r41_22-pic0-11.json': '03ff00ffaaa2f197276a6b1018e45f8cf3fc3a19a1371e8352ebb6c3fabeeae5', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r41_22-pic0-12.json': '0841229ace4fb69381af167994ea78c01e6cfec128921accab7d536255d05aac', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r41_22-pic0-13.json': 'd9bf37bf1aea9d9371f7875edbf01b019a72349054a80c98f35444ad5f9869a0', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r41_22-pic0-14.json': '1ebf15c2bd5115665a7303ca649c39a1b95435f6db78a51b957a22b72cf8a2de', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r41_22-pic0-15.json': 'b9e64cf25fd87fe31f280223f337a8efecf44363e3f0210e70c8add0e5baddeb', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r41_22-pic0-32.json': '457f6e063becc706d72f28da45f45e8ef95ced1dc239b9b660be302dfaa52679', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r41_22-pic0-33.json': '321f3d26f12697f6a846a4fe485523afd919889374b591a8cb8988c2f507e7e2', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r41_22-pic0-34.json': '001e4065c95c4241fb8e18f5580bc810945be79c6e8bf98ad5ac6f1d6c765896', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r41_22-pic0-35.json': 'ffa4867e334ed4eec3a2e9d0c1abcd7b855dfb2140848cffbcae91977e9c5eae', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r41_22-pic0-36.json': '1fe60ba4e206032a41c87c7f9e57de84b83e5ab90ed3095fb8f3f8512274104f', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r41_22-pic0-37.json': '4c6b23ee60b1b7615588d3f26f064fc4fb41bc5d6ee377aa90fa1fbf59933d26', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r41_22-pic0-38.json': 'e9f3ef81961fc679aedc88b9a5da95319c3d89b0abff1eee69dc289d4f79b1dd', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r41_22-pic0-39.json': '101bbf631c12a55e7ad667b47ced4db14c01bdd7ee00188a73304117f481b352', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r41_22-pic0-40.json': '74b87389edaa927caeb03a95162eaabaf03cf51c4c1a85ac9a50662112d960b1', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r41_22-pic0-41.json': 'd9fa5514b7ba5b1c4040117dfb124408a746a7d70f010a33db6dc81972ff92ad', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r41_22-pic0-42.json': '46962e0f991fd151adf0a4fe853c5e4ea16d4a44d29f797444421cbe93d80458', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r41_22-pic0-43.json': '0f7c84706c9a0f28fe17f18a38132ff4d4d9ba962eaf3d263684e59511af623d', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r41_22-pic0-44.json': 'fb88903c7908c2536e639b91d2505c7fbc11c61be8878250b4706f1c75631fb6', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r41_22-pic0-45.json': '3d3a14f7f426978e57827a0fe2c1f9c676d54195b0e014805eba50d6f2345083', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r41_22-pic0-46.json': '16fa0b75f0aa228c8ede07b48969e4b369ae426f2a403733b66f81cb774266e9', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r41_22-pic0-47.json': 'b4ac450a798e026524a994941e419b7b239315075fbb6271130c6b3e7d795e9f', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r43_20-pic0-00.json': 'e0ed71b4a8f4c3d4fd46efcec72512f05dd9e6547e0a23a35bc7dddb45cc269f', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r43_20-pic0-01.json': 'f23ba2cf1ffbd6dfdbe45a53062ee4508875e3008cfc2cc50e765b6030e62278', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r43_20-pic0-02.json': '91d05e9f380f9de6243806a3d1b1062f00d763783074c33b101bb8cd38fa1f23', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r43_20-pic0-03.json': 'f0595776ab9d0e6750c4b8d0045c5c149421123e67c56629ea37f42282c8dfe9', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r43_20-pic0-04.json': '0b3d341cd694e55ab1c953b92f21b7a73143db0dc6278cb98cea855e3bda9be2', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r43_20-pic0-05.json': 'fd07bc484d26f93d68449f3c8d65610acdb37954d406be5795376c5bf50363a2', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r43_20-pic0-06.json': '57fddc07cd0fe0141afbc32dd01e76f85d0ad0a06770dc6f14eb5cc4f8f299b4', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r43_20-pic0-07.json': 'b7364f0cb1d7fe579ac259ab07c73b5e2112cc4c8d7eac01091200d0724a1e59', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r43_20-pic0-08.json': 'd194b63bfd576ce9ffc65829c7514a30ce4d85360dce4ecfd0d60d9b67a5bdb5', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r43_20-pic0-09.json': 'd482a24d144a7dadb32af3f7cbdc43c25abf8d75ed35c68c3fa527a02b1ad7f1', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r43_20-pic0-10.json': '1931de9e72e5293f0f3ee1ba3919814666ffaccb834453d2b8064101b0e6e624', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r43_20-pic0-11.json': '1b7af7b29fafc890fef5224289db2d71d3ea3ec033c952041dfe74dc43b875f0', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r43_20-pic0-12.json': 'da8bebed9080c5c3c49b1044448c282b89846b001a731cd67e783b1f70ca03d1', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r43_20-pic0-13.json': '25debf81ee7fccef76977c642799e302f79ba42dfbac916bb21e70d900270d58', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r43_20-pic0-14.json': '5a84f4dcbf294da86f2db53128bc053eedcb754e35e8fe93f1d3ad5c1f7af0ae', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r43_20-pic0-15.json': 'e398aaf0e1d26a4aae8385de9fbba1da028099ee52b1006596bd58c87ac3fefd', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r43_20-pic0-32.json': '6b8640bf651783ede7c0191e12162f02c5a64e624af6bf4f2416961b7d44602b', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r43_20-pic0-33.json': '3027fe5c783a9af06981948b8325bc098caf7a2286768697ce4c9b53132c5106', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r43_20-pic0-34.json': 'f7e50b4684163b77c78d0b53bc90a63c0de3837d13280a735cc7bc56ff815b6e', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r43_20-pic0-35.json': 'ca4265d2ae11e299a33f0c6f66208d75c032d05f0b5f459bd824257ea2e43b18', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r43_20-pic0-36.json': '02f686bb5870828291bea84a1ab600136d00985b27f86e1154b7b7d6ecf1afa2', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r43_20-pic0-37.json': '6be5e62bca1ea4f43b610145f7cc5d25163f9643ea3daa1e865472dbfd3fd13b', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r43_20-pic0-38.json': '798b3d986825708e007265ef9a1a17c7e7c456a77e8c7e03d3a72427050380f5', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r43_20-pic0-39.json': 'a688f890fc3fbd53610115130d93010ea2303dad38301d1b45c8fd68769b3f82', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r43_20-pic0-40.json': 'd05db6284c7b6a621568d0da344b0996d19b8774d2fd744286d4cbfc163b8670', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r43_20-pic0-41.json': '0212bb58408e60be764646d045bf3192c271a359de200eb926e0592c6730c11f', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r43_20-pic0-42.json': 'cf840dd115f9f755d1294eeda238f5777820d99c11f3620a0b8e509a69033a3a', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r43_20-pic0-43.json': '01a51a30c97e9348d3fbcbc62cb2130f533a74a28f192a1c18602e1ae8f84629', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r43_20-pic0-44.json': 'b0fa20462c4e66d9dd7ef8f15f88ae9a2bedd8885c5e9030a149f9da15a4d017', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r43_20-pic0-45.json': 'd4b5eefe035035fb2429d1d79f7acf3ad55ea6bfeb401b5d98d5af6e0d912293', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r43_20-pic0-46.json': 'f3f2b2fab099d483d7078c8bb269cb440a2b9cbe4a374e2e67f04aff4d00c455', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r43_20-pic0-47.json': '7158e46cb020e7aed6532c143886f6593b3cea780c8c5ad38552ce2fbfdc36b8', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r43_31-pic0-00.json': 'd5761d20e40b535b87dab4763a0c22c72b88d8d706759bc8aa04b5ca80e6dacb', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r43_31-pic0-01.json': '3575f0011ba4b5fe436989a330f5a20a06193956e3daf032d2a097dbf0d9f8e2', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r43_31-pic0-02.json': '1aca7d5713a7e8c00ac6d17774605d71199fbf160449dc4d953f0d54ac48eb60', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r43_31-pic0-03.json': 'f551701fc602bd678790fdddc598223102ac6ecd3d71e7ef1dd55a52ec694db4', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r43_31-pic0-04.json': '725640fd8d032cd21075853cfc6d392e409dac54048049122219eaaf9ef0bdc4', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r43_31-pic0-05.json': '524a85b035b2a93d624703229a2037a8409079cceedb1824cf71af0b21b8aab7', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r43_31-pic0-06.json': 'b30474610a4d961ac6f949e0b8bc56ad55561da4980d1ff4bad54a6e7cba130b', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r43_31-pic0-07.json': '4e7b0ae4c09087e4191dcc4eaefa41fdcd4048f8fd410fecfca158efd56b77e2', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r43_31-pic0-08.json': '8e655e812c66bef9eec450be73859297ab8c0e283c0c465f66000624bdb340ef', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r43_31-pic0-09.json': 'e24261290695af6ba1f7d3acc2ffff93962b3bb8ee6df50b1c2076f0658e8459', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r43_31-pic0-10.json': '992fa66d7e64646d8f1d43b93cba25efc5aa36c3f906dd932586538a8b2afc5f', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r43_31-pic0-11.json': '415156563f73d701ae5ebbb966df93b26ecafd25e2e4c9f460fd057fbfbe861e', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r43_31-pic0-12.json': '6a1afd43699a6e6f51556d242357e8c6c5650f6e22b5e5ea11a5f4c859aaabe0', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r43_31-pic0-13.json': '866dc7c584b949398d3db3a4e838b25f7810587e999b0fa4b2d40de27131bba6', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r43_31-pic0-14.json': 'f9fc4ea79895a9cceeca85ec7179cf3ca895b76b8d64dc05df1a70318573abcb', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r43_31-pic0-15.json': '864e959ceaf53df828ff83291fb7a318743ef6110fef364c372040e40c15c24d', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r43_31-pic0-32.json': '207926dc5d1fd55a1d2f125a667b91424018ccc1d7fa394293219e9c8b637a00', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r43_31-pic0-33.json': 'b25352dd42cadc1947a77afd4deae8fa6a719e2b6b7f02521433390a0ac1f68c', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r43_31-pic0-34.json': '668d6af8cf0b4f4a2696eaa47424f9d9af384d0e4456e89363c789638d9eb620', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r43_31-pic0-35.json': '429acd870dc12a2452ff6f763873177477216ff06016ed88a3bee9c34f9ae884', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r43_31-pic0-36.json': '30f928c6b1462c0653d9c012830e23b3a74588e628a780a419e259448505a731', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r43_31-pic0-37.json': 'f6f6d8be1f0c68e7d0ad56297a4d4cdedad76cf97ed5836da22e565f3151b291', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r43_31-pic0-38.json': 'b316718e4c7c22929bdea838c1b7f87f41da88942a6b58c0d66107667e594c2c', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r43_31-pic0-39.json': '03136aa577be9966c05cef441f5be07a8c18fbeb058f1e23798527483dab20cf', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r43_31-pic0-40.json': '21d1a82f3bef6e96ed8e88bf64ed581e2be4024e5ed52e8750d6841a55c83af7', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r43_31-pic0-41.json': '2ef8460ffc27f4df233f386f8ed155dbf0a587a02a307d81a193dbfc2cd35d69', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r43_31-pic0-42.json': 'e33af35f886ebe0d55c116970b14a77be095380aab741cab8d9885105ffbca38', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r43_31-pic0-43.json': 'ca28417ad0f80c868d2fc4bff69c47965f88a8af195c147a58e07fb2b0248ada', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r43_31-pic0-44.json': 'b8232f6addd935804927532e675db084d65a7f4798ccbf013b618d36bc1c06fe', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r43_31-pic0-45.json': '88bc45ec20619c03d5b5f5a2e933d91e360e5f6625ff3e7a2a85907d51c7a7e4', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r43_31-pic0-46.json': '8aa0a1574b9c355f6213d67eacf35d1ad90f3cb0ea20c5631bc93050cdfd5a88', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r43_31-pic0-47.json': 'ff0c5bf5728d5f714d5ccf669c8fc0e65b3be751004094fcf5954f539a61bef4', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r7-pic0-00.json': '71ee6367726f7a88f40c86321365249b2ec36f193b612c21813c39b5820eb848', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r7-pic0-01.json': '7b9919dd544d8dfcff65ddb449d1601d5b00b49051cbe15e6d7af317f084ff7f', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r7-pic0-02.json': '1c6fd85291bf6eac46e9daff0ea0a12a1038e6057c4770619c8f14d93303f728', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r7-pic0-03.json': '6b34134af61ad20a0b3024d6f91edd1cc5426f5e1bd664ffdf0a1edfc449e672', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r7-pic0-04.json': '72717797afad1b76b314c155e943f5b90677c6e87d791b7da13c25a06e9fb077', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r7-pic0-05.json': '10cbc7fc7aaf2bd11a5d58f06291540e449b2a656305dcc8c1e2b1aea8e3070b', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r7-pic0-06.json': '585cdc9458e0b00cd4f592e34e4004a382d72ed1ebbca571e3378f6b4ab03c74', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r7-pic0-07.json': '6ef7c6ef8644bfabe8608e9a52ec4b0c4b624659c954f05a9575dc5cf481a443', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r7-pic0-08.json': '47b8f13200effa824f26bf0e93c73b0d584a6d7bd683f245e308860c7a2b2f82', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r7-pic0-09.json': '945c97aab3368f745d9b56b7bdf31f67244ed1521bfc91f7a421e3f27256c173', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r7-pic0-10.json': '4f072343ad90c8ce5ed134b6e472307fb546448b1497d0ef57468e99e95d38d2', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r7-pic0-11.json': '5cb8a984e40c017b5e278b07338da7014e4eb4348393cf377e56418ec1cd23a3', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r7-pic0-12.json': 'fabb3623fe3c51aba4cc4f5e2d2437f5f6665035dd3ac7124fe740376d195d75', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r7-pic0-13.json': '5b3a6297cc57f7715c9129343c68d7583ed15070bf7a6909504226f14d320f2e', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r7-pic0-14.json': '8be9d3f368e9aa2386bdcb86ad84c22ca5820fce82a0671a3fe1ddc005daffe2', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r7-pic0-15.json': 'f60123e30e248fe07406c4df5ab1ac9972a716027d9e5af80917dd46d4568ff6', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r7-pic0-32.json': 'c38b64d858291c20cccc14cbe4ebfef6494162b980b34ab19eb65c818854251e', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r7-pic0-33.json': '3d571cb031446813c9576a765ffd438330bece12f319a023e7c85ecdb4690edf', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r7-pic0-34.json': 'aa6d7ea324bf6ff421c4786f3f4305a28c084ce0dc47c915050c9dec1a553d4c', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r7-pic0-35.json': '670603504ecebad45e9f5371a997a83f5e5854520abfefab92949505dba51388', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r7-pic0-36.json': 'b6fb8bf480312f8a1a1f4d5ead5fb7e69cd64bbde7cdbc8a18c6341d897d7758', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r7-pic0-37.json': '573cd74d9a2b189ea57acbf110f7bceeae55bb7903a5ba873b98bd55af7b9948', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r7-pic0-38.json': '2e61f1d834950083bc2015d892649e22199115266e1dbe4dfb8e97e38d7bc449', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r7-pic0-39.json': '064a8924b40901ce2b5afb6735a441382d30753af723e4f7a0f3cc37265bc5de', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r7-pic0-40.json': 'ccb82ca904d8a9d5cc7b129193d9f382a628a53a1894dc2dbc11f3f28d3a2723', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r7-pic0-41.json': '3c894567e76a67e96fb2c682cab0f2dcdc5184b8456af45f23872ab46fb02402', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r7-pic0-42.json': 'c3f08da22c1bb315b2437710eb5104e13984393023f4185299aceae770723e14', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r7-pic0-43.json': '8bb8d8fe9f40f858da705208145e3f510e96d74ecbd948ecc6679c36999beb30', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r7-pic0-44.json': '2a9426cf13c745b30d30266740296251bcecc591a9f5c12747ed3b39829042ae', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r7-pic0-45.json': '8c64f0a5a0f8c9a1db7dc81c2f65e28a3e226f9700173f0ff493cd9ef48cf1f7', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r7-pic0-46.json': '3833cf16a1f255c101ee8fda025f538a9386d79e15bfcc988e3ba08a03671da3', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r7-pic0-47.json': '8deba1c2d59843ae88f3bc6145ab5c18fdf8d84f611eebf86207c88469078c05', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r7sq-pic0-00.json': '587c19f0042ca30a6478624eedb5f781ebee237dcba4807871a0de0c984e5348', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r7sq-pic0-01.json': '5c9c09554a26ca0560affdac11cebd7a8dd0005d09d7ea324dea30419d2d94ca', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r7sq-pic0-02.json': '755ebd305b9e63b38302eada234a75f5a4de240272949c5489f99fc4a9ecd535', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r7sq-pic0-03.json': 'b378b64eaf123ad2f9c3a7e604f0b0b1093eaea651ac22a851ae3cf2e8cd07c7', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r7sq-pic0-04.json': 'd7b91d9f26fa909c9276e79f2e25e37865079d263728d8b3026ea8ab0ae82b4b', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r7sq-pic0-05.json': 'c3af6539d0815f56dddd4cf4b8879340207b41f218f85401d7d9eefcedac1e7a', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r7sq-pic0-06.json': '90fbd369c88001cd06e6436438670576c8757e9fc8339eabf0fa793bd082153e', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r7sq-pic0-07.json': '3a1af18e3a7b3e70666c8f8545a07ea4af0dd8ddb5af34fc3b6a45a19c351414', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r7sq-pic0-08.json': '5c5e05f5fd9446ac1a4a1f82343d7d8218d48834a3113576c40b185ed63092f0', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r7sq-pic0-09.json': 'cd2b5383896122d844ed2a6d22d6fea38fcdef7eb732f3a9eb29297954b05d39', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r7sq-pic0-10.json': 'd43bbfc105209e59df3072b4e361bb793abb38ed9ca470b76b2c5e3c7210e8fa', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r7sq-pic0-11.json': '7d2a00d1e22d230973aa8dcf908e0fff2fe8eb7aeb9fd90b3051574c8745a694', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r7sq-pic0-12.json': '5adc22e17b0c91a7d096c225b03b9e2d77dc0602846fea4a3633a3e080dbf375', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r7sq-pic0-13.json': '84d867abdbe78562720eb2f2e838126a16e88bb7b78e5f6b1936409e9a54e4a1', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r7sq-pic0-14.json': 'f34834d314c167ff35e4a8629240cbd61501070478feb0a8abb21e76f657f530', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r7sq-pic0-15.json': 'b5d7d9491bf48e39fcfcf8fda49ee14c42db0abd662bce935075460a14497d64', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r7sq-pic0-32.json': '194725c4575508cab61e9eb00161e28f34daf53d609c9c10dd02fb7d25309304', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r7sq-pic0-33.json': 'e8990aaabbe3c28fc7485e3df7d310b7e1a0afbeb0f4cd824df65fabdd51f28d', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r7sq-pic0-34.json': 'e77c968730cc6d69d07d6ce77c5a8f196123ed017bd52aa6fd7f349d9168d90a', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r7sq-pic0-35.json': '516102cc3e4066026ea48d12b8a3e369854c8eef9b035e2dced6e066747d746f', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r7sq-pic0-36.json': '4087f926f0f006dc0eda96e40b36a426f8689530d83b55f8e341fc57dd291115', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r7sq-pic0-37.json': '9b735e46bf1164521b1ecaeabc4a6b413f2a0f54edeb2322b594d9874dad4001', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r7sq-pic0-38.json': '1000715b526b0065d8c0f4135191fa57d2ce43769d16e5226ab0c158a63c5235', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r7sq-pic0-39.json': '914a2d8a0b87122abf811e7f44af7d67d8fa048b32fb8232aa79af44b436ea07', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r7sq-pic0-40.json': '6b235b6cbce1eca4df03bdb6e580ae811103a03d667607d5b2f2a70e9d946479', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r7sq-pic0-41.json': '6c8e4747e2b041aead67041a1bb195c6a184dabd19c3c6ee06ddb2c9c5ae9c07', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r7sq-pic0-42.json': '5866d537fa2b4417b554041a9783b629154bc0e8545252a347b88e9f73ba25c3', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r7sq-pic0-43.json': '6b69d656555b55858d3f55d0e0e1f161f4298ccd7d75580e0414d48b24b5c43c', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r7sq-pic0-44.json': '9d6583d589fb3e72f6215c2193713db5e23444eaaa36dcefe6bebbb9f67d94de', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r7sq-pic0-45.json': 'd961658d9cd01cbf30b334b206cfb95a454ca795fee02256b93d593a4511e6bd', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r7sq-pic0-46.json': '54c62d6d3680ccd02706a37e0320b752435eb73da9da234fc385d75a933e90cd', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r7sq-pic0-47.json': '2862cf24c52eca213437619702595922afa6f819fbc0954d44ebc208d3bcdd29', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r83-pic0-00.json': '9661a334505a7ffaa0390f5eaa455dba0a2ab5298be42bc4cab60e9b5a5dec0c', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r83-pic0-01.json': 'ccbf703792a44487dd187379952237d2b5bbe86391a8f4d01bae39350e7072ea', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r83-pic0-02.json': '0b1ea9263bb7cfd984fa3130798ca7afa034e1ff245218ae798a9bb11aada5fc', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r83-pic0-03.json': '56cca01870a9ff8fa7c6d53dc4a56c965c24686f06192518186a8568b82e22b4', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r83-pic0-04.json': '1a0ea57a4feea80c7fd4d6d709a53d3f83dadd3348658079bdb417702cca1a9a', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r83-pic0-05.json': '55cad79b756fd9a4d449b1a92bc4b3444aa3c6227a83582d5451b8c6b439866c', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r83-pic0-06.json': '055ef1d84333f1d14b65948166e487379374639024c8fa8506ed1492756a62a4', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r83-pic0-07.json': '9e229492efcf970d02afb33295d8fa3e8875a84c0c193ba0e53cd73e464d7025', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r83-pic0-08.json': 'a57b4629c5dd09b7543de354a1c0ed624d3b51a55a3b32f7f555958fb08d679a', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r83-pic0-09.json': 'c386912ab3187d2981613c652e41b71ac2539c1396ebd439a0115d6f80b9f737', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r83-pic0-10.json': '07a680ab2e5e75da80ec23da23b04ad5ed35e10fd19653785505426c07d3420b', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r83-pic0-11.json': '78fb94840ac13a376b9e34970b74eb3df84126d8aebbe10a8b221b80dfab31fe', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r83-pic0-12.json': '91bf8ec21f643545d7b59633c69881272f8fde71fcfc1cd83140980a2a0de107', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r83-pic0-13.json': '775be3ea5041ac768599d431a9f8c3f2d88d25ab162c24bae0236d908844d5e7', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r83-pic0-14.json': '93da8c84e9ce680630ecac9d1237b5191619bc667556e5c56910db433d3d3b74', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r83-pic0-15.json': 'c3b0f937fba4e74273cc4c9eabf8a6bd155e54d39dde9f9a2f430a7c5aa44637', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r83-pic0-32.json': 'bf9880173a1b657d6bfcaa42fff041ba3e4820588890af95d6ef10efb27eece5', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r83-pic0-33.json': 'a1807121dc8c515a05a6adcf21260d380cc572f719576e7be242f5623324afa2', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r83-pic0-34.json': 'db586ce05a20ba1bc711847a140e5de12924c85e8f43d5d7f27addeb344f2e9c', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r83-pic0-35.json': 'aada25fb3da8cfa30ae973ebb578b21993897b57ee8b2bc1dd0ac389925998c9', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r83-pic0-36.json': 'add79fe69f1d078872ab81d8922ab2ed098994356403767376043ae976192e9d', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r83-pic0-37.json': '6f9ac56300157c9a8c2ef8b598f1201942d306f821f08d8739267c40e1f0e9bb', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r83-pic0-38.json': 'e0263c3ec7c89b7393a59efeaf81bc44cc3368076c97c2d1d66e9bac20fa5083', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r83-pic0-39.json': '9db559522040f17666632b2937b7b3fcae7c36b14c87c186410f2c9d389e8de2', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r83-pic0-40.json': 'e7c33b9f586e305356fd88ba76920745bc3b9491ad435af8c0221a3a017aa94d', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r83-pic0-41.json': 'f41b1311107d5fa4ecc6a9ed029ff3f15fe8b7dc8c16f0b5804cd07ee2cf149a', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r83-pic0-42.json': '1d0c81e07a9e7764eb32d037cef13fef636bd8fed93c9905cf1e7a56e9849603', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r83-pic0-43.json': '39ec705592c0f007ae31a09774cde9176b58cb8980c62a396d813c632022d677', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r83-pic0-44.json': '9c6e19dfa73e1c6310cba5682cfdec5f1c70ff5b3836cfbd602dbba0836355e5', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r83-pic0-45.json': 'c2a91f26a70f057123e84ae41cffe94b2da0676f01c0f001a90d79304d51c2fb', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r83-pic0-46.json': '298bc434d6b05356713425d736479ce1ebf830f0a473247d764232dc5eeb4248', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r83-pic0-47.json': '0a3cc5c01272b1a7cf195110d138be1d7a82d14872d8588da1cf256adff6697f', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r83-pic1-00.json': 'e645491baaaf0cdedb6ea040f0ee7ee7b010e50fd5ae25270a5732809361e3e7', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r83-pic1-01.json': '7110936233a31bed08baf9e9ab692002232f69d9ff42ffcb0e1615d7f0e33403', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r83-pic1-02.json': '5e366cb6dcec8ca3c2e09fe1770e0e8feb4c83fe513e0c47a7fafe571b98f683', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r83-pic1-03.json': 'a200f0d0b8da8e41665f27016d711c39ed5875663ae75a2834bf8e5f15141078', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r83-pic1-04.json': '004307002b6776a45184c90f71425762b47a58ff3a6009e9c6d8dcf07ef05808', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r83-pic1-05.json': '04d3b3d95c3d3c367296abf0fad96fa6b26ede546849b834017ee37a53a46255', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r83-pic1-06.json': '4970a439cc00ac8d219ade18834b8541deff7cdcb314b9788d3184a0f0785841', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r83-pic1-07.json': 'da98b87cf9e93b1d62a2f4cbbc622025d0e404a84f8f9a79dbd173da1d59a162', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r83-pic1-08.json': 'f5aee264d6154ab2e3c9d350d3634d11ad840a3d0bf681eef1615a0c5f419bea', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r83-pic1-09.json': 'c4c7cd783b965452afe7335458de0531ef03b2e62cb76e5c3305dc2b00fe1d23', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r83-pic1-10.json': '63a447f19f3d1c10031b63e3000b6b0a7e11937da9b10f0c5c2ad8dd2750d302', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r83-pic1-11.json': '02d5db073a5c64af3ed65d669c9a7d0a81fd49c38031384a25e34538fe2f3dd5', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r83-pic1-12.json': 'ccee2ea1d8b19c88fffb09fc7cc87d9c1d0f5bb19849b9d42bd622168da47615', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r83-pic1-13.json': 'c126e0cf71b15d5edab2de7899c2e59d7281e0f1e9cad0103c2bb8a70fdc00a2', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r83-pic1-14.json': '2cbfdbaddc72265bfb6bbf1f0d7a33697e2195fbd4c1b2480acf7facce97a9b6', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r83-pic1-15.json': 'e35976107ab0d213164deb9402a755fd0a1dc5c396cd9abda52517bdd854286b', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r83-pic1-16.json': 'a782f63e0d8c3b2981cc402a2203a7b540aaa11a30a5d2ead616b79fb2504184', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r83-pic1-17.json': '3dc796228b136be4e0a7e58208f4030bd06fcb4dc38f15f33ee7799841f36de6', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r83-pic1-18.json': '6238f40262649a8403eb11629cd559b0fb53230185620fd9b66624a51e43ac5f', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r83-pic1-19.json': 'f49d40d1cb57a422a6e132b4e80a9eb8db56941679a33071ed5bcba997eae835', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r83-pic1-20.json': '8287f64739b0a4614540e1eb2d62300851b2d5ef01e1b53e0d8ad16c8de61aad', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r83-pic1-21.json': '4a0ccfa5878eff27e66c96e5fef9cbbd3c0fe80ef187bd1f721b124c76d6ae9f', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r83-pic1-22.json': 'aa9f59ce7df191815a36013bdd674df00a99ab95ccbb308716549a07bf8b5c9f', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r83-pic1-23.json': '0dce7d945a705dd0ab5ae9526e3afdc65335cbbe78f5f136421f9386795d5540', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r83-pic1-24.json': '8d3bcb27aa813aefc2b334c566bbb27fa7ed98070a2cd9c9da1376b61918612c', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r83-pic1-25.json': 'c03111cd61bec8e7251205938beb034b3c2c65fe269f5d08c71e11e0ba29725c', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r83-pic1-26.json': '62df6ffc28929d286104b4239480631066a38489289bd99c2e5d2243dd039594', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r83-pic1-27.json': 'c5576ba30c8b552cea230496c46fd8442195184ba14e27427ab1a1f5cf8bf909', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r83-pic1-28.json': 'dc6e898550ed6d5215ba15e2d711f1e9c132592a3061cf852a50d70b08bf6fe6', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r83-pic1-29.json': 'f3ed95e9c79800df1c8ad507ae421f21c2a9babb229083d08e39cdb8f3c66528', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r83-pic1-30.json': '8226d8e1b4cff257e99aa182001907318474394cf7d878abd2e47276901e8ac0', 'certificates/reflected_packet_v83/plans/reflected_ab_v1/r83-pic1-31.json': 'f6722b707b3c15c03e465df7f032912aa1e1022c5f38c5fd0777e25c2acc1bda', 'certificates/reflected_packet_v83/positive-witness-review.json': '2d512c0e16d1aebb10ec5fe543b4243aae8ee926d14800ef186345e3d71aea63', 'certificates/reflected_packet_v83/positive-witnesses.json': '3698b92f9dd7bcc6176d4845515e1723263548fe6ad4a9d93451e15972707518', 'certificates/reflected_packet_v83/primary_target.py': 'b94aca2bc29fc612fbecaadf83430ebd0a9e09f462876c1dc41fe0aa51422bd3', 'certificates/reflected_packet_v83/selection.json': 'b9514f79d11a0e1e520d2ff2f49dfaca8ef0ddf72bf4c91c98fa2b8b0d6dbc65', 'certificates/reflected_packet_v83/source-relocation.json': '969f0ca3022f15fc351ce161e5930e87d8b92f89849e4d47596a020654d65bda', 'certificates/reflected_primary_kernel_v83.py': 'f76493968959d467ff661d7a54cbadaa51fdf0b9e9fee2fbd2968160193d69b6', 'certificates/residual_orbit_k17_k38_v42.py': '5186066f21a2f840ffc71936a54372d7ad51eb6cc919bf7812bb1bec258cd59d', 'independent/d7_ambient_points_sturm_v45.py': '3eb68f57e6d0d88b2e548f605742d434c8714b25eb1404b3c27d32415b61f42a', 'independent/d7_r83_principal_sturm_v55.py': '894339e82389ff651bf25b04693821c4dfeedc427dec2fa4fd8993f195bb663f', 'independent/quadrilateral_k17_stdlib_v41.py': 'c5bd45f334a028820a071b5c758acfb8e48afc1e8454b18a3174dd04ca58a38c', 'independent/reflected_packet_v83/sturm_target.py': '6a0aae7a310a091719721b08654574da87f0c105e716988e932ceedacd5886e6', 'independent/reflected_sturm_kernel_v83.py': 'b2b4b154198290c66f7f1c41b639a23bb62943a64772146c2cd11ca6cb5a1c4e', 'independent/residual_orbit_k17_k38_sturm_v42.py': '12a0d875e2715fe2c2ae123d2dc5202f4d7ce4e024f114ad163126654219fc13'}
DATA = 'certificates/reflected_packet_v83'
RELOCATION_SHA256 = '969f0ca3022f15fc351ce161e5930e87d8b92f89849e4d47596a020654d65bda'
MATH_PINS = {'certificates/d7_ambient_points_v45.py': '34cb35db91b953162b595fb61e031f551e6df0a0d7b5192dc15a6c00c39274a2', 'certificates/d7_r83_mask14_v46.py': '70e9c7114adea2b319c3e430046ae3e917717142e7a78e4d43a043950f1147f3', 'certificates/d7_r83_masks14_30_v47.py': 'e78ba5e793134f28d8113d8129809421a5a7754c439afcd64ba0385543799941', 'certificates/d7_r83_principal_expected_v55.py': '0fd36a67b004e6a250188aa8731a207456b5e49faa955005d6b12cd36665bea7', 'certificates/d7_r83_principal_v55.py': 'a937b69037b5d53f40a600aefb177bf89766d3529b22f73529b2b861ec9238cc', 'certificates/full_packet_v62/ambient_points.json': 'c0e84da6c266c087cbce65a6e8aca8e3f30350ce56b5a1089c4ae4978a381f36', 'certificates/full_packet_v62/module_templates.json': '09ccff217511a67e946567135f9ffb74da020843621d0895ecce8746d25c245d', 'certificates/reflected_packet_v83/primary_target.py': 'b94aca2bc29fc612fbecaadf83430ebd0a9e09f462876c1dc41fe0aa51422bd3', 'certificates/reflected_primary_kernel_v83.py': 'f76493968959d467ff661d7a54cbadaa51fdf0b9e9fee2fbd2968160193d69b6', 'certificates/residual_orbit_k17_k38_v42.py': '5186066f21a2f840ffc71936a54372d7ad51eb6cc919bf7812bb1bec258cd59d', 'independent/d7_ambient_points_sturm_v45.py': '3eb68f57e6d0d88b2e548f605742d434c8714b25eb1404b3c27d32415b61f42a', 'independent/d7_r83_principal_sturm_v55.py': '894339e82389ff651bf25b04693821c4dfeedc427dec2fa4fd8993f195bb663f', 'independent/quadrilateral_k17_stdlib_v41.py': 'c5bd45f334a028820a071b5c758acfb8e48afc1e8454b18a3174dd04ca58a38c', 'independent/reflected_packet_v83/sturm_target.py': '6a0aae7a310a091719721b08654574da87f0c105e716988e932ceedacd5886e6', 'independent/reflected_sturm_kernel_v83.py': 'b2b4b154198290c66f7f1c41b639a23bb62943a64772146c2cd11ca6cb5a1c4e', 'independent/residual_orbit_k17_k38_sturm_v42.py': '12a0d875e2715fe2c2ae123d2dc5202f4d7ce4e024f114ad163126654219fc13'}
EXECUTED_CAMPAIGN_PINS = {'certificates/d7_ambient_points_v45.py': '34cb35db91b953162b595fb61e031f551e6df0a0d7b5192dc15a6c00c39274a2', 'certificates/d7_r83_mask14_v46.py': '70e9c7114adea2b319c3e430046ae3e917717142e7a78e4d43a043950f1147f3', 'certificates/d7_r83_masks14_30_v47.py': 'e78ba5e793134f28d8113d8129809421a5a7754c439afcd64ba0385543799941', 'certificates/d7_r83_principal_expected_v55.py': '0fd36a67b004e6a250188aa8731a207456b5e49faa955005d6b12cd36665bea7', 'certificates/d7_r83_principal_v55.py': 'a937b69037b5d53f40a600aefb177bf89766d3529b22f73529b2b861ec9238cc', 'certificates/full_packet_v62/ambient_points.json': 'c0e84da6c266c087cbce65a6e8aca8e3f30350ce56b5a1089c4ae4978a381f36', 'certificates/full_packet_v62/module_templates.json': '09ccff217511a67e946567135f9ffb74da020843621d0895ecce8746d25c245d', 'certificates/residual_orbit_k17_k38_v42.py': '5186066f21a2f840ffc71936a54372d7ad51eb6cc919bf7812bb1bec258cd59d', 'experiments/continuous_v83/child.py': 'e593f070efa845a05d9f7d4ab7d1efd2566ebc275bad8a9e02621a56f0548360', 'experiments/continuous_v83/nominate.py': '2b91049e3e37609d4c1a94edd230c8f2fb5730a9ab49ad4c9018b5aa87fe8f89', 'experiments/continuous_v83/primary.py': 'c167a51e6881923319576897ea07a8d05998f935bc0f6884eef529928c4d6313', 'experiments/continuous_v83/primary_target.py': 'b94aca2bc29fc612fbecaadf83430ebd0a9e09f462876c1dc41fe0aa51422bd3', 'experiments/continuous_v83/runner.py': '3acda0f59005ad92ae832760e4c614c524eed078e84b571698f8209ff5614449', 'experiments/continuous_v83/sturm.py': '9d83a1bcaa79ebbbd4fed557cd7845de57e3955e6f55516ee7119a425f45d2ca', 'experiments/continuous_v83/sturm_target.py': '6a0aae7a310a091719721b08654574da87f0c105e716988e932ceedacd5886e6', 'independent/d7_ambient_points_sturm_v45.py': '3eb68f57e6d0d88b2e548f605742d434c8714b25eb1404b3c27d32415b61f42a', 'independent/d7_r83_principal_sturm_v55.py': '894339e82389ff651bf25b04693821c4dfeedc427dec2fa4fd8993f195bb663f', 'independent/quadrilateral_k17_stdlib_v41.py': 'c5bd45f334a028820a071b5c758acfb8e48afc1e8454b18a3174dd04ca58a38c', 'independent/residual_orbit_k17_k38_sturm_v42.py': '12a0d875e2715fe2c2ae123d2dc5202f4d7ce4e024f114ad163126654219fc13'}
TARGETS = ('ORIGINAL_Q_V1', 'REFLECTED_AB_V1')
POSITIVES = {'OM:pic0:07': 'A', 'OM:pic0:15': 'B', 'OM:pic0:00': 'C'}
CANONICAL = ('preparation.json', 'A-entries.json', 'C-entries.json', 'bottom-rows.json',
             'candidates-nonzero.json', 'candidates-zero.json', 'classification.json', 'survivors.json')
COUNTS = ('a_count', 'c_count', 'bottom_rows', 'nonzero_candidates', 'zero_candidates',
          'candidate_count', 'survivor_count', 'classification_digest', 'c_zero_audit')
ZERO = (0,) * 6


def require(value, message):
    if not value:
        raise ValueError(message)


def canonical(value):
    return (json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False) + '\n').encode()


def object_hash(value):
    return hashlib.sha256(canonical(value)).hexdigest()


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def safe(path):
    path = Path(os.path.abspath(path))
    for component in (path, *path.parents):
        require(not component.is_symlink(), 'Symlink path: ' + str(component))
    return path


def at(root, relative):
    require(type(relative) is str and relative and not relative.startswith('/') and '\\' not in relative
            and all(part not in ('', '.', '..') for part in relative.split('/')), 'Unsafe relative path')
    return safe(root / relative)


def regular(path):
    info = safe(path).stat()
    require(stat.S_ISREG(info.st_mode) and info.st_nlink == 1, 'Nonregular/aliased evidence file: ' + str(path))
    return info


def no_duplicates(items):
    result = {}
    for key, value in items:
        require(key not in result, 'Duplicate JSON key: ' + key)
        result[key] = value
    return result


def read(path):
    regular(path)
    with Path(path).open() as stream:
        return json.load(stream, object_pairs_hook=no_duplicates,
            parse_constant=lambda value: (_ for _ in ()).throw(ValueError('Nonfinite JSON')))


def write(path, value):
    path = safe(path)
    with path.open('xb') as stream:
        stream.write(canonical(value)); stream.flush(); os.fsync(stream.fileno())
    return digest(path)


def valid_hash(value):
    return type(value) is str and re.fullmatch('[0-9a-f]{64}', value) is not None


def checked_sources():
    require(BACKEND in ('primary', 'sturm') and KERNEL ==
            ('certificates/reflected_primary_kernel_v83.py' if BACKEND == 'primary' else 'independent/reflected_sturm_kernel_v83.py'),
            'Frontend backend/kernel contract')
    require(len(MATH_PINS) == 16 and len(EXECUTED_CAMPAIGN_PINS) == 19, 'Source provenance sizes')
    require(all(SOURCE_PINS.get(name) == checksum for name, checksum in MATH_PINS.items()),
            'Missing or changed relocated mathematical closure')
    for name in ('selection.json', 'expected.json', 'source-relocation.json'):
        require(DATA + '/' + name in SOURCE_PINS, 'Missing frozen data pin: ' + name)
    for relative, expected in SOURCE_PINS.items():
        require(relative.startswith(('certificates/', 'independent/')) and valid_hash(expected), 'Nonportable source/data pin')
        path = at(ROOT, relative); regular(path)
        require(digest(path) == expected, 'Source/data changed: ' + relative)
    require(SOURCE_PINS[DATA + '/source-relocation.json'] == RELOCATION_SHA256, 'Unreviewed relocation manifest')
    return dict(SOURCE_PINS)


def source_mapping():
    record = read(at(ROOT, DATA + '/source-relocation.json'))
    require(type(record.get('schema_version')) is int and record['schema_version'] == 1 and
            record['portable_sha256'] == MATH_PINS and
            record['executed_campaign_source_sha256'] == EXECUTED_CAMPAIGN_PINS and
            record['mathematical_file_count'] == 16, 'Original19/portable16 source mapping mismatch')
    # The old v55 expected module is an unchanged, pinned import-lineage file.
    # It is never parsed here as this run's expected outcome.
    return record


def literal_pins(path):
    values = [node.value for node in ast.parse(path.read_text()).body if isinstance(node, ast.Assign)
              and any(isinstance(t, ast.Name) and t.id == 'SOURCE_PINS' for t in node.targets)]
    require(len(values) == 1, 'Kernel source pin declaration')
    return ast.literal_eval(values[0])


def declared_labels():
    rows = read(at(ROOT, 'certificates/full_packet_v62/ambient_points.json'))
    labels = {row['point_label'] for row in rows}
    require(len(rows) == len(labels) == 368 and set(POSITIVES) <= labels, 'C39 declaration changed')
    return labels


def validate_plan(plan, target, label):
    require(type(plan) is dict and set(plan) == {'schema_version', 'target_id', 'point_label', 'transport', 'entry_units'}, 'Plan schema')
    require(type(plan['schema_version']) is int and plan['schema_version'] == 2 and
            plan['target_id'] == target and plan['point_label'] == label, 'Wrong plan target/label')
    route = plan['transport']; entries = plan['entry_units']
    require(type(route) is dict and set(route) == {'inversion_first', 'unit_exponents'} and
            type(route['inversion_first']) is bool and type(entries) is dict and set(entries) == {'A', 'C', 'DEN', 'NUM'},
            'Plan route/unit schema')
    for powers in [route['unit_exponents'], *entries.values()]:
        require(type(powers) is list and len(powers) == 5 and
                all(type(x) is int and -2 <= x <= 2 for x in powers), 'Unit nomination exponent domain')


def select_job(target, label):
    selection = read(at(ROOT, DATA + '/selection.json'))
    require(type(selection.get('schema_version')) is int and selection['schema_version'] == 1 and
            selection.get('scope') == 'PORTABLE_FULL_TARGET_REPLAY_SELECTION_NO_PROMOTION', 'Selection schema/scope')
    require(type(selection.get('jobs')) is dict and selection['jobs'], 'Missing exact replay job domain')
    labels = declared_labels()
    for key, row in selection['jobs'].items():
        require(type(row) is dict and set(row) == {'target_id', 'point_label', 'kind', 'plan_path', 'plan_sha256', 'cpu_seconds'},
                'Selection job schema')
        require(row['target_id'] in TARGETS and row['point_label'] in labels and
                key == row['target_id'] + '/' + row['point_label'], 'Selection identity')
        require(row['kind'] == ('KNOWN_POSITIVE_CONTROL' if row['point_label'] in POSITIVES else 'NEGATIVE_REPLAY'), 'Job kind/label mismatch')
        require(row['plan_path'].startswith(DATA + '/plans/') and
                SOURCE_PINS.get(row['plan_path']) == row['plan_sha256'], 'Frozen plan path/seal missing')
        require(type(row['cpu_seconds']) is int and 1 <= row['cpu_seconds'] <= 7200, 'Unbounded CPU selection')
    key = target + '/' + label
    require(target in TARGETS and key in selection['jobs'], 'Target/label outside exact frozen replay selection')
    return selection, selection['jobs'][key]


def installed_runtime():
    """Hash package sources without importing SymPy into the independent process."""
    result = {'invocation_path': sys.executable, 'resolved_binary': str(Path(sys.executable).resolve()),
              'python_sha256': digest(Path(sys.executable).resolve()), 'python_version': sys.version,
              'prefix': sys.prefix, 'base_prefix': sys.base_prefix, 'platform': platform.platform(),
              'machine': platform.machine(), 'packages': {}}
    for name, version in (('sympy', '1.14.0'), ('mpmath', '1.3.0')):
        require(importlib.metadata.version(name) == version, 'Unreviewed package version: ' + name)
        spec = importlib.util.find_spec(name)
        require(spec is not None and spec.origin is not None, 'Missing installed package: ' + name)
        directory = Path(spec.origin).resolve().parent
        files = {p.relative_to(directory).as_posix(): digest(p) for p in sorted(directory.rglob('*'))
                 if p.is_file() and p.suffix in ('.py', '.so', '.pyd', '.dylib')}
        require(files, 'Empty runtime source inventory')
        result['packages'][name] = {'version': version, 'files': files}
    return result


def cpu_limit(job):
    wanted = job['cpu_seconds']; limits = resource.getrlimit(resource.RLIMIT_CPU)
    if limits == (resource.RLIM_INFINITY, resource.RLIM_INFINITY):
        resource.setrlimit(resource.RLIMIT_CPU, (wanted, wanted + 1))
        limits = resource.getrlimit(resource.RLIMIT_CPU)
    require(0 < limits[0] <= wanted and limits[0] <= limits[1] <= wanted + 1,
            'Inherited CPU limits exceed the frozen replay cap')
    return list(limits)


def invoke_kernel(target, label, output, plan):
    spec = importlib.util.spec_from_file_location('portable_full_v83_' + BACKEND, at(ROOT, KERNEL))
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    old_bytecode = sys.dont_write_bytecode; sys.dont_write_bytecode = True
    try:
        spec.loader.exec_module(module)
        if BACKEND == 'primary':
            return module.run_point(label, output, False, plan, target_id=target)
        old_argv = sys.argv
        sys.argv = [str(at(ROOT, KERNEL)), '--target-id', target, '--point', label,
                    '--output', str(output), '--witness-plan', str(plan)]
        try:
            module.main()
        finally:
            sys.argv = old_argv
    finally:
        sys.dont_write_bytecode = old_bytecode


def output_inventory(directory):
    inventory = {}; inodes = set()
    for parent, children, files in os.walk(safe(directory), followlinks=False):
        for name in children:
            safe(Path(parent) / name)
        for name in files:
            path = safe(Path(parent) / name); info = regular(path)
            require((info.st_dev, info.st_ino) not in inodes, 'Duplicate output inode')
            inodes.add((info.st_dev, info.st_ino))
            inventory[path.relative_to(directory).as_posix()] = {'size_bytes': info.st_size, 'sha256': digest(path)}
    return inventory


def flat(value, width):
    require(type(value) is list and len(value) == width, 'Bad canonical matrix/row width')
    answer = []
    for vector in value:
        require(type(vector) is list and len(vector) == 6 and
                all(type(x) is str and re.fullmatch('0|-?[1-9][0-9]*', x) for x in vector), 'Bad integral coefficient vector')
        answer.extend(map(int, vector))
    return tuple(answer)


def ordered(values, width):
    answer = set(); previous = None
    require(type(values) is list, 'Canonical list required')
    for value in values:
        item = flat([value] if width == 1 else value, width)
        require(previous is None or previous < item, 'Nonunique/noncanonical numeric order')
        answer.add(item); previous = item
    return answer


def sign_closed(values):
    require({tuple(-x for x in item) for item in values} == values, 'Incomplete sign restoration')


def inspect_artifacts(output, result):
    a = ordered(read(output / 'A-entries.json'), 1); c = ordered(read(output / 'C-entries.json'), 1)
    bottoms = ordered(read(output / 'bottom-rows.json'), 2)
    nz = ordered(read(output / 'candidates-nonzero.json'), 4); z = ordered(read(output / 'candidates-zero.json'), 4)
    for rows in (a, c, bottoms, nz, z):
        sign_closed(rows)
    require(ZERO in a and ZERO in c and ZERO + ZERO in bottoms and
            all(row[:6] in c for row in bottoms), 'Original domain membership/origins')
    require(all(g[12:18] != ZERO for g in nz) and all(g[12:18] == ZERO for g in z) and nz.isdisjoint(z), 'Wrong determinant branch')
    candidates = nz | z
    require(all(g[:6] in a and g[12:] in bottoms for g in candidates), 'Candidates outside original domains')
    classification = read(output / 'classification.json')
    require(set(classification) == {'excluded', 'survivors', 'histogram'}, 'Classification schema')
    survivors = ordered(read(output / 'survivors.json'), 4); sign_closed(survivors)
    require(ordered(classification['survivors'], 4) == survivors, 'Survivor lists differ')
    excluded = set(); histogram = Counter(); digest_value = hashlib.sha256(); previous = None
    for row in classification['excluded']:
        require(set(row) == {'gamma', 'failure'}, 'Classification row schema')
        gamma = flat(row['gamma'], 4); failure = row['failure']
        require(previous is None or previous < gamma, 'Classification order/duplicate')
        require(type(failure) is list and len(failure) == 2 and type(failure[0]) is int and
                0 <= failure[0] < 4 and type(failure[1]) is int and 0 <= failure[1] < 6, 'Polygon index')
        excluded.add(gamma); previous = gamma; histogram[tuple(failure)] += 1
        digest_value.update((','.join(map(str, gamma)) + ';' + ','.join(map(str, failure)) + '\n').encode())
    require(excluded.isdisjoint(survivors) and excluded | survivors == candidates, 'Incomplete classification')
    require(classification['histogram'] == [{'side_index': s, 'embedding_index': e, 'count': n}
            for (s, e), n in sorted(histogram.items())], 'Classification histogram')
    audit = result['c_zero_audit']
    require(type(audit) is list and len(audit) == 5 and all(type(v) is int and v >= 0 for v in audit), 'c=0 audit shape')
    require(audit[0] == sum(row[:6] == ZERO for row in bottoms) and audit[1] == 1 and
            audit[1] + audit[2] <= audit[0] and audit[3] <= audit[2] and audit[4] == len(z), 'c=0 audit counts')
    counts = dict(zip(COUNTS, (len(a), len(c), len(bottoms), len(nz), len(z), len(candidates),
                              len(survivors), digest_value.hexdigest(), audit)))
    require(all(result.get(k) == v and (type(result[k]) is int if type(v) is int else True)
                for k, v in counts.items()), 'Counts disagree with complete canonical artifacts')
    checkpoint = set()
    with (output / 'survivors-partial.jsonl').open() as stream:
        for line in stream:
            row = json.loads(line, object_pairs_hook=no_duplicates)
            require(set(row) == {'gamma', 'status'} and row['status'] == 'INCOMPLETE_STAGE_SURVIVOR_CHECKPOINT', 'Survivor checkpoint schema')
            gamma = flat(row['gamma'], 4); require(gamma not in checkpoint, 'Repeated survivor checkpoint'); checkpoint.add(gamma)
    require(checkpoint == survivors, 'Survivor checkpoints incomplete')
    return counts, candidates, survivors


def positive_witnesses(key, expected, plan_hash, counts, candidates, survivors):
    paths = [DATA + '/positive-witnesses.json', DATA + '/positive-witness-review.json']
    require(all(path in SOURCE_PINS for path in paths), 'Missing reviewed positive witness inputs')
    data, review = [read(at(ROOT, path)) for path in paths]
    require(type(data.get('schema_version')) is int and data['schema_version'] == 1 and
            data.get('status') == 'REVIEWED_EXACT_POSITIVE_WITNESSES', 'Positive witness data status')
    require(type(review.get('schema_version')) is int and review['schema_version'] == 1 and
            review.get('status') == 'PASS_EXACT_POSITIVE_WITNESS_REVIEW' and
            review.get('source_relocation_sha256') == RELOCATION_SHA256, 'Positive witness review status/source')
    row = data['records'][key]
    require(type(row) is dict and {'target_id', 'point_label', 'plan_sha256', 'preparation_sha256',
            'target_geometry_sha256', 'signed_H_selected', 'evidence_sha256'} <= set(row) and
            type(row['evidence_sha256']) is dict and row['evidence_sha256'] and
            all(valid_hash(value) for value in row['evidence_sha256'].values()), 'Missing positive witness evidence bindings')
    checksum = object_hash(row)
    require(review['record_sha256'].get(key) == checksum == expected.get('positive_witness_record_sha256'),
            'Positive witness record lacks exact reviewed binding')
    require(row['target_id'] == expected['target_id'] and row['point_label'] == expected['point_label'] and
            row['plan_sha256'] == plan_hash and row['target_geometry_sha256'] == expected['target_geometry_sha256'] and
            row['preparation_sha256'] == expected['canonical_artifact_sha256']['preparation.json'], 'Witness selected-coordinate identity')
    signed = ordered(row['signed_H_selected'], 4)
    require(len(signed) == 2 and {tuple(-x for x in g) for g in signed} == signed and
            signed <= candidates and signed <= survivors, 'Reviewed signed H missing from candidates or survivors')
    require(counts['survivor_count'] >= 2, 'Known positive did not retain its signed witnesses')
    return {'record_sha256': checksum, 'review_sha256': SOURCE_PINS[paths[1]], 'both_signed_witnesses_present': True}


def validate_complete(result, output, expected, job, plan_hash, runtime, limits, inventory):
    target, label = job['target_id'], job['point_label']; key = target + '/' + label
    require(result.get('status') == 'FULL_GEOMETRY_DISCOVERY_COMPLETE' and result.get('target_id') == target and
            result.get('point_label') == label and result.get('backend') == BACKEND, 'Fresh full target/label result required')
    require(all(result.get(k) is True for k in ('assertions_enabled', 'sources_unchanged', 'executable_unchanged', 'witness_plan_unchanged')),
            'Missing kernel guards')
    require(result.get('driver_sha256') == MATH_PINS[KERNEL] and
            result.get('source_sha256') == literal_pins(at(ROOT, KERNEL)), 'Wrong kernel source closure')
    raw_runtime = result.get('runtime', {})
    require(raw_runtime.get('python_executable') == runtime['resolved_binary'] and
            raw_runtime.get('python_sha256') == runtime['python_sha256'] and
            raw_runtime.get('python_version') == runtime['python_version'] and
            raw_runtime.get('cpu_limits_seconds') == limits, 'Kernel runtime/cap identity mismatch')
    for name in ('cpu_seconds', 'wall_seconds'):
        value = result.get(name)
        require(type(value) in (int, float) and math.isfinite(value) and value >= 0, 'Invalid kernel timing')
    require(result.get('input_witness_plan_sha256') == plan_hash == result.get('witness_plan_sha256'), 'Kernel did not consume the exact frozen plan bytes')
    observed = {name: inventory[name]['sha256'] for name in CANONICAL}
    require(observed == result.get('canonical_artifact_sha256'), 'Kernel canonical file inventory mismatch')
    artifacts = result.get('artifact_sha256')
    require(type(artifacts) is dict and artifacts, 'Missing declared artifact seals')
    for name, checksum in artifacts.items():
        at(output, name)
        require(valid_hash(checksum) and inventory.get(name, {}).get('sha256') == checksum, 'Declared output seal mismatch: ' + name)
    require(all(artifacts.get(name) == checksum for name, checksum in observed.items()), 'Declared artifact closure omits canonical data')
    target_hash = inventory['target-geometry.json']['sha256']
    require(result.get('target_geometry_sha256') == target_hash == artifacts.get('target-geometry.json') and
            read(output / 'target-geometry.json').get('target_id') == target, 'Target descriptor mismatch')
    require(expected['target_id'] == target and expected['point_label'] == label and expected['kind'] == job['kind'] and
            expected['plan_sha256'] == plan_hash and expected['target_geometry_sha256'] == target_hash,
            'Frozen expected target/point/plan mismatch')
    require(observed == expected['canonical_artifact_sha256'], 'Eight canonical files differ from reviewed expectations')
    require(type(expected.get('counts')) is dict and set(expected['counts']) == set(COUNTS), 'Expected nine-count schema')
    require(all(type(expected['counts'][name]) is int and expected['counts'][name] >= 0 for name in COUNTS[:7]) and
            valid_hash(expected['counts']['classification_digest']) and
            type(expected['counts']['c_zero_audit']) is list and len(expected['counts']['c_zero_audit']) == 5 and
            all(type(v) is int and v >= 0 for v in expected['counts']['c_zero_audit']), 'Malformed expected counts')
    counts, candidates, survivors = inspect_artifacts(output, result)
    require(counts == expected['counts'], 'Nine complete counts differ from reviewed expectations')
    require(type(expected.get('evidence_sha256')) is dict and expected['evidence_sha256'] and
            all(valid_hash(value) for value in expected['evidence_sha256'].values()), 'Missing frozen outcome provenance')
    positive = None
    if job['kind'] == 'NEGATIVE_REPLAY':
        require(not survivors and counts['survivor_count'] == 0, 'Survivor cannot pass an exclusion replay')
    else:
        positive = positive_witnesses(key, expected, plan_hash, counts, candidates, survivors)
    return observed, counts, positive


def execute(target, label, output):
    require(__debug__, 'Assertions must be enabled')
    output = safe(output)
    require(not output.exists(), 'Refusing to overwrite any output directory')
    started = time.monotonic(); own_hash = digest(Path(__file__))
    sources = checked_sources(); source_mapping()
    _, job = select_job(target, label)
    plan_path = at(ROOT, job['plan_path']); plan_hash = digest(plan_path)
    plan = read(plan_path); validate_plan(plan, target, label)
    limits = cpu_limit(job)
    runtime_before = installed_runtime()
    try:
        # No new expected geometry or positive-witness record is parsed before
        # this unchanged kernel completes ALL enumeration/classification stages.
        invoke_kernel(target, label, output, plan_path)
        result = read(output / 'result.json')
        require(result.get('status') == 'FULL_GEOMETRY_DISCOVERY_COMPLETE', 'Nonfull kernel result cannot be a mathematical pass')
        expected_table = read(at(ROOT, DATA + '/expected.json'))
        require(type(expected_table.get('schema_version')) is int and expected_table['schema_version'] == 1 and
                expected_table.get('scope') == 'REVIEWED_COMPLETE_TARGET_OUTCOMES', 'Expected table schema/scope')
        expected = expected_table['points'][target + '/' + label]
        inventory = output_inventory(output)
        observed, counts, positive = validate_complete(result, output, expected, job, plan_hash, runtime_before, limits, inventory)
        require(read(output / 'witness-plan.json') == plan, 'Emitted plan differs from frozen target-aware data')
        runtime_after = installed_runtime()
        require(runtime_after == runtime_before, 'Runtime/interpreter/package sources changed during full replay')
        require(checked_sources() == sources and digest(Path(__file__)) == own_hash and digest(plan_path) == plan_hash,
                'Frontend/source/plan changed during full replay')
        require(output_inventory(output) == inventory, 'Kernel output changed during frontend validation')
        inventory_hash = write(output / 'kernel-output-inventory.json', inventory)
        before_hash = write(output / 'frontend-runtime-before.json', runtime_before)
        after_hash = write(output / 'frontend-runtime-after.json', runtime_after)
        record = {'schema_version': 1, 'status': 'POINT_PASS', 'target_id': target, 'point_label': label,
            'backend': BACKEND, 'kind': job['kind'], 'scope': 'FRESH_FULL_TARGET_REPLAY_NO_AUTOMATIC_PROMOTION',
            'identity_key': target + '/' + label + '/' + BACKEND,
            'frontend_sha256': own_hash, 'source_sha256': sources,
            'portable_math_source_sha256': MATH_PINS, 'executed_campaign_source_sha256': EXECUTED_CAMPAIGN_PINS,
            'source_relocation_sha256': RELOCATION_SHA256, 'kernel_result_sha256': inventory['result.json']['sha256'],
            'kernel_output_inventory_sha256': inventory_hash, 'canonical_artifact_sha256': observed, 'counts': counts,
            'target_geometry_sha256': inventory['target-geometry.json']['sha256'],
            'input_plan_sha256': plan_hash, 'emitted_plan_sha256': inventory['witness-plan.json']['sha256'],
            'expected_record_sha256': object_hash(expected), 'expected_file_sha256': SOURCE_PINS[DATA + '/expected.json'],
            'runtime_before_sha256': before_hash, 'runtime_after_sha256': after_hash,
            'approved_runtime_acceptance': 'EXTERNAL_SUPERVISOR_REVIEW_REQUIRED',
            'invocation_path': sys.executable, 'resolved_binary': runtime_before['resolved_binary'],
            'cpu_limits_seconds': limits, 'sources_unchanged': True, 'runtime_unchanged': True,
            'executable_unchanged': True, 'witness_plan_unchanged': True, 'assertions_enabled': __debug__,
            'positive_witness_check': positive, 'wall_seconds_including_frontend': time.monotonic() - started}
        write(output / 'certificate-result.json', record)
        print('POINT_PASS ' + json.dumps(record, sort_keys=True), flush=True)
        return record
    except BaseException as error:
        # Never erase, replace or relabel kernel outputs after a failed run.
        if not output.exists():
            output.mkdir(parents=True)
        failure = {'schema_version': 1, 'status': 'FULL_REPLAY_FAILED_NO_MATHEMATICAL_PASS',
                   'target_id': target, 'point_label': label, 'backend': BACKEND,
                   'error_type': type(error).__name__, 'error': str(error), 'frontend_sha256': own_hash,
                   'input_plan_sha256': plan_hash, 'wall_seconds': time.monotonic() - started}
        write(output / 'frontend-failure.json', failure)
        raise


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--target-id', choices=TARGETS, required=True)
    parser.add_argument('--point', required=True)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args(argv)
    # Retain temporary outputs even on failure; a suite must preserve their path.
    output = args.output if args.output is not None else Path(tempfile.mkdtemp(prefix='portable-full-v83-')) / 'artifacts'
    print('RETAINED_ARTIFACT_DIRECTORY ' + str(safe(output)), flush=True)
    return execute(args.target_id, args.point, output)


if __name__ == '__main__':
    main()
